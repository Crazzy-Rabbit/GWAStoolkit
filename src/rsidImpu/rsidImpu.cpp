//
//  rsidImpu.cpp
//  GWAStoolkit
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/linereader.hpp"
#include "utils/writer.hpp"  // support .gz format out
#include "utils/log.hpp"
#include "utils/util.hpp"
#include "utils/gwasQC.hpp" // basic QC
#include "utils/FormatEngine.hpp"
#include "rsidImpu/rsidImpu.hpp"
#include "rsidImpu/allele.hpp"
#include "rsidImpu/dbsnp.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <deque>
#include <vector>

using namespace std;

// 结构体：保存排序后的 GWAS 信息（不改变原始行顺序）
struct GWASRecord {
    size_t index;         // 原始在 gwas_lines 中的下标
    string chr;           // 规范化后的染色体
    long long pos;        // 位置
    string a1;            // 原始 A1
    string a2;            // 原始 A2
    string canon_a1;      // 规范化后的 A1
    string canon_a2;      // 规范化后的 A2
};

void process_rsidImpu(const Args_RsidImpu& P)
{

    deque<string> gwas_lines;
    string line;
    
    //================ 1. 读取 GWAS header =================
    LineReader reader(P.gwas_file);
    if (!reader.getline(line)) {
        LOG_ERROR("Empty GWAS summary file.");
        exit(1);
    }

    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);

    // chech exist of SNP col
    bool has_SNP = false;
    int idx_SNP  = -1;

    for (size_t i = 0; i < header.size(); i++){
        string hlow = header[i];
        std::transform(hlow.begin(), hlow.end(), hlow.begin(), ::tolower);

        if (hlow == "snp"){
            has_SNP = true;
            idx_SNP = i;
            break;
        }
    }
    
    // find the col, CHR, POS, A1, A2, P 列
    int gCHR     = find_col(header, P.g_chr); 
    require(gCHR >= 0, "GWAS missing required column [" + P.g_chr + "] for rsidImpu.");
    
    int gPOS     = find_col(header, P.g_pos);
    require(gPOS >= 0, "GWAS missing required column [" + P.g_pos + "] for rsidImpu.");

    int gA1      = find_col(header, P.g_A1);
    require(gA1 >= 0, "GWAS missing required column [" + P.g_A1 + "] for rsidImpu.");

    int gA2 = find_col(header, P.g_A2);
    require(gA2 >= 0, "GWAS missing required column [" + P.g_A2 + "] for rsidImpu.");

    int idx_beta = find_col(header, P.col_beta);
    int idx_se   = find_col(header, P.col_se);
    int idx_freq = find_col(header, P.col_freq);
    int idx_pv   = find_col(header, P.g_p);
    int idx_n    = find_col(header, P.col_n);

    //================ 2. 读入 GWAS 数据行 =================
    while (reader.getline(line)){
        if (line.empty()) continue;
        gwas_lines.push_back(line);
    }
    
    size_t n = gwas_lines.size();
    LOG_INFO("Loaded GWAS lines (data): " + std::to_string(n));

    //================ 基础 QC：过滤无效 N/beta/se/freq/P =================
    vector<bool> keep_qc(n, true);
    double maf   = P.maf_threshold;

    bool can_qc = (idx_beta >= 0 &&
                idx_se   >= 0 &&
                idx_freq >= 0 &&
                idx_pv   >= 0 &&
                idx_n    >= 0);
    if (can_qc) {
        gwas_basic_qc(gwas_lines, header, idx_beta, idx_se, idx_freq, idx_pv, idx_n, keep_qc, maf);
    } else {
        LOG_WARN("Cannot perform full QC in rsidImpu (missing beta/se/freq/N/p columns).");
    }

    //================ 构建排序用的 GWASRecord 向量 =================
    vector<GWASRecord> gwas_vec;
    gwas_vec.reserve(n);

    for (size_t i = 0; i < n; ++i){
        auto f = split(gwas_lines[i]);

        GWASRecord rec;
        rec.index = i;
        rec.chr   = canonical_chr(f[gCHR]);
        rec.pos   = std::stoll(f[gPOS]);
        rec.a1    = f[gA1];
        rec.a2    = f[gA2];

        auto canon   = canonical_alleles(rec.a1, rec.a2);
        rec.canon_a1 = canon.first;
        rec.canon_a2 = canon.second;

        gwas_vec.push_back(std::move(rec));
    }

    // 按 chr, pos 排序（稳定排序）
    std::stable_sort(gwas_vec.begin(), gwas_vec.end(),
        [](const GWASRecord& a, const GWASRecord& b){
            if (a.chr != b.chr) return a.chr < b.chr;
            return a.pos < b.pos;
        }
    );

    LOG_INFO("GWAS records sorted by CHR:POS for two-pointer matching.");

    //================ 准备匹配结果容器 =================
    vector<bool>   keep(n, false);   // 是否最终进入主输出
    vector<string> rsid_vec(n);      // 匹配到的 rsID

    //================ 单通扫描 dbSNP（Two-pointer merge） =================
    LineReader dbr(P.dbsnp_file);
    string dline;

    bool is_bim = ends_with(P.dbsnp_file, ".bim") ||
                    ends_with(P.dbsnp_file, ".bim.gz");

    int dCHR, dPOS, dA1, dA2, dRS;
    vector<string> dhdr;

    if (!is_bim) {
        // 有 header 的一般表格格式
        if (!dbr.getline(dline)) {
            LOG_ERROR("Empty dbSNP file.");
            exit(1);
        }
        dhdr = split(dline);

        dCHR = find_col(dhdr, P.d_chr);
        dPOS = find_col(dhdr, P.d_pos);
        dA1  = find_col(dhdr, P.d_A1);
        dA2  = find_col(dhdr, P.d_A2);
        dRS  = find_col(dhdr, P.d_rsid);

        if (dCHR<0 || dPOS<0 || dA1<0 || dA2<0 || dRS<0){
            LOG_ERROR("dbSNP header incomplete.");
            exit(1);
        }
    } else {
        // .bim / .bim.gz 格式：CHR RSID CM POS A1 A2
        dCHR = 0; dRS = 1; dPOS = 3; dA1 = 4; dA2 = 5;
    }

    LOG_INFO("Start two-pointer merge between GWAS and dbSNP.");

    size_t gi = 0;
    const size_t Gn = gwas_vec.size();
    size_t scanned = 0;

    while (dbr.getline(dline)){
        if (dline.empty()) continue;
        auto f = split(dline);

        string dchr = canonical_chr(f[gCHR]);
        long long dpos = std::stoll(f[dPOS]);

        // 推进 GWAS 指针：直到 gwas_chr:pos >= dbsnp_chr:pos
        while (gi < Gn && 
                (gwas_vec[gi].chr < dchr ||
                (gwas_vec[gi].chr == dchr && gwas_vec[gi].pos < dpos))) {
            ++gi;
        }

        if (gi >= Gn) break;

        // 如果 chr 不同，继续读 dbSNP
        if (gwas_vec[gi].chr > dchr) {
            continue;
        }

        // 现在有两种可能：
        // 1) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos == dpos → 候选匹配
        // 2) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos > dpos → 说明 dbSNP 这行在 GWAS 中不存在该 pos，继续读下一行 dbSNP
        if (!(gwas_vec[gi].chr == dchr && gwas_vec[gi].pos == dpos)) {
            continue;
        }

        // 等位基因规范化
        auto canon_db = canonical_alleles(f[dA1], f[dA2]);
        const string& db_a1 = canon_db.first;
        const string& db_a2 = canon_db.second;
        const string& rsid  = f[dRS];

        // 可能有多个 GWAS 行在同一 chr:pos（或者多个 dbSNP 行同一 chr:pos）
        // 对所有该位置的 GWAS 进行尝试匹配
        size_t gj = gi;
        while (gj < Gn &&
                gwas_vec[gj].chr == dchr &&
                gwas_vec[gj].pos == dpos) {

            size_t orig_idx = gwas_vec[gj].index;

            // QC 未通过的行不做 rsID 匹配，但仍保留为 keep=false（之后输出到 unmatch）
            if (keep_qc[orig_idx]) {
                if (gwas_vec[gj].canon_a1 == db_a1 &&
                    gwas_vec[gj].canon_a2 == db_a2) {
                    // 找到匹配
                    keep[orig_idx]     = true;
                    rsid_vec[orig_idx] = rsid;
                }
            }
            ++gj;
        }

        ++scanned;
        if(scanned % 1000000 == 0){
            LOG_INFO("[dbSNP two-pointer] scanned " + std::to_string(scanned/1000000) + "M lines.");
        }
    }

    LOG_INFO("Two-pointer merge finished. dbSNP lines scanned: " + std::to_string(scanned));

    //================ 去重（按 rsID / P 值） =================
    if (P.remove_dup_snp) {
        gwas_remove_dup(
            gwas_lines,
            header,
            idx_pv,          // P 列 index
            rsid_vec,
            keep
        );
    }
    
    //================ Writer（自动 txt / gz） =================
    bool out_is_gz = ends_with(P.out_file, ".gz");

    std::string out_main    = P.out_file;
    std::string base        = P.out_file;
    if (out_is_gz && base.size() > 3) {
        base.erase(base.size() - 3);  // 去掉 ".gz"
    }
    std::string out_unmatch = out_is_gz ?
                                (base + ".unmatch.gz") :
                                (P.out_file + ".unmatch");
    
    Writer fout(out_main, P.format);
    Writer funm(out_unmatch, P.format);

    if (!fout.good() || !funm.good()) {
        LOG_ERROR("Error opening output file.");
        exit(1);
    }
    
    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);

    // header
    if (P.format == "gwas") {
        // 原始 header + SNP
        string h;
        for (size_t j=0; j<header.size(); j++) {
            if (j) h += "\t";
            h += header[j];
        }
        if (!has_SNP){
            h += "\tSNP";
        }
        fout.write_line(h);
    } else {
        // 格式化 header
        string h;
        for (size_t j=0; j<spec.cols.size(); j++) {
            if (j) h += "\t";
            h += spec.cols[j];
        }
        fout.write_line(h);
    }

    // outfile 
    for (size_t i=0; i<n; i++){
        gwas_lines[i].erase(
            std::remove(gwas_lines[i].begin(), gwas_lines[i].end(), '\r'),
            gwas_lines[i].end()
        );

        // 未保留（QC 不通过或没匹配到 rsID）→ unmatch
        if (!keep[i]) {
            funm.write_line(gwas_lines[i]);
            continue;
        }

        auto f = split(gwas_lines[i]);

        if (P.format == "gwas") {
            if (has_SNP){
                // ★ 覆盖原 N
                if (idx_SNP < (int)f.size()){
                    f[idx_SNP] = rsid_vec[i];
                }
                // 重建一行
                string out;
                for (size_t j=0; j<f.size(); j++){
                    if (j) out +="\t";
                    out += f[j];
                }
                fout.write_line(out);
            } else {
                // ★ 原来没有 SNP → 追加
                string out = gwas_lines[i] + "\t" + rsid_vec[i];
                fout.write_line(out);
            }
            
            continue;
        }

        // 构造 row 映射 (col→value)
        unordered_map<string,string> row;
        row["SNP"]  = rsid_vec[i];
        row["A1"]   = f[gA1];
        row["A2"]   = f[gA2];

        if (idx_freq>=0 && idx_freq<(int)f.size())
            row["freq"]=f[idx_freq];
        if (idx_beta>=0 && idx_beta<(int)f.size())
            row["beta"]=f[idx_beta];
        if (idx_se>=0 && idx_se<(int)f.size())
            row["se"]=f[idx_se];
        if (idx_pv>=0 && idx_pv<(int)f.size())
            row["p"]=f[idx_pv];
        if (idx_n>=0 && idx_n<(int)f.size())
            row["N"]=f[idx_n];

        fout.write_line(FE.format_line(spec, row));
    }
}