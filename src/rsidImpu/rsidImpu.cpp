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

#include <fstream>
#include <iostream>
#include <algorithm>
#include <deque>
#include <vector>
#include <unordered_map>

using namespace std;

// =======================================================
// GWAS minimal record for merge
// =======================================================

struct GWASRecord {
    size_t    index;     // original line index
    int       chr;       // chr_code
    int64_t   pos;
    AlleleKey allele;
};

// 用 tab 定位，第 col_idx 列替换为 value
static inline void replace_nth_column_inplace(
    std::string &line,
    int col_idx,
    const std::string &value
){
    size_t start = 0;
    for (int c = 0; c < col_idx; ++c){
        start = line.find('\t', start);
        if (start == std::string::npos) return;
        ++start;
    }

    size_t end = line.find('\t', start);
    if (end == std::string::npos){
        // 最后一列
        line.replace(start, line.size() - start, value);
    } else {
        line.replace(start, end - start, value);
    }
}

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

    bool can_qc = (idx_beta >= 0 ||
                idx_se   >= 0 ||
                idx_freq >= 0 ||
                idx_pv   >= 0 ||
                idx_n    >= 0);
    if (can_qc) {
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(gwas_lines, header, idx_beta, idx_se, idx_freq, idx_pv, idx_n, keep_qc, maf);
    } else {
        LOG_WARN("Cannot perform full QC in rsidImpu (missing beta/se/freq/N/p columns).");
    }

    //================ 构建排序用的 GWASRecord 向量 =================
    vector<GWASRecord> gwas_vec;
    gwas_vec.reserve(n);

    for (size_t i = 0; i < n; ++i){
        auto f = split(gwas_lines[i]);

        int chr = canonical_chr_code(f[gCHR]);
        if (chr < 0) continue;

        int64_t pos = std::strtoll(f[gPOS].c_str(), nullptr, 10);
        if (pos <= 0) continue;

        AlleleKey ak = make_allele_key(f[gA1], f[gA2]);
        if (ak.type == 2) continue;

        GWASRecord rec;
        rec.index  = i;
        rec.chr    = chr;
        rec.pos    = pos;
        rec.allele = ak;

        gwas_vec.push_back(rec);
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

    size_t gi = 0, Gn = gwas_vec.size();
    size_t scanned = 0;

    while (dbr.getline(dline)){
        if (dline.empty()) continue;
        auto f = split(dline);

        int dchr = canonical_chr_code(f[dCHR]);
        if (dchr < 0) continue;

        int64_t dpos = std::strtoll(f[dPOS].c_str(), nullptr, 10);

        // 推进 GWAS 指针：直到 gwas_chr:pos >= dbsnp_chr:pos
        while (gi < Gn &&
            (gwas_vec[gi].chr < dchr ||
            (gwas_vec[gi].chr == dchr && gwas_vec[gi].pos < dpos))) {
            ++gi;
        }

        if (gi >= Gn) break;
        // 现在有两种可能：
        // 1) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos == dpos → 候选匹配
        // 2) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos > dpos → 说明 dbSNP 这行在 GWAS 中不存在该 pos，继续读下一行 dbSNP
        if (gwas_vec[gi].chr != dchr || gwas_vec[gi].pos != dpos) continue;

        // 等位基因规范化
        AlleleKey db_allele = make_allele_key(f[dA1], f[dA2]);
        if (db_allele.type == 2) continue;

        const string& rsid = f[dRS];
        // 可能有多个 GWAS 行在同一 chr:pos（或者多个 dbSNP 行同一 chr:pos）
        // 对所有该位置的 GWAS 进行尝试匹配
        size_t gj = gi;
        while (gj < Gn && 
            gwas_vec[gj].chr == dchr &&
            gwas_vec[gj].pos == dpos) {

            size_t orig_idx = gwas_vec[gj].index;
            // QC 未通过的行不做 rsID 匹配，但仍保留为 keep=false（之后输出到 unmatch）
            if (keep_qc[orig_idx] &&
                gwas_vec[gj].allele.type == db_allele.type &&
                gwas_vec[gj].allele.key  == db_allele.key) {

                // 正向匹配 || 反向匹配
                keep[orig_idx]     = true;
                rsid_vec[orig_idx] = rsid;
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

        if (P.format == "gwas") {
            if (has_SNP){
                // ★ 覆盖原 N
                replace_nth_column_inplace(
                    gwas_lines[i],
                    idx_SNP,
                    rsid_vec[i]
                );
                fout.write_line(gwas_lines[i]);
            } else {
                // 原本没有 SNP：直接追加
                fout.write_line(gwas_lines[i] + "\t" + rsid_vec[i]);
            }

            continue;
        }

        // 构造 row 映射 (col→value)
        auto f = split(gwas_lines[i]);

        unordered_map<string,string> row;
        row["SNP"]  = rsid_vec[i];
        row["A1"]   = f[gA1];
        row["A2"]   = f[gA2];

        if (idx_freq>=0 && idx_freq<(int)f.size()) row["freq"]=f[idx_freq];
        if (idx_beta>=0 && idx_beta<(int)f.size()) row["beta"]=f[idx_beta];
        if (idx_se>=0 && idx_se<(int)f.size()) row["se"]=f[idx_se];
        if (idx_pv>=0 && idx_pv<(int)f.size()) row["p"]=f[idx_pv];
        if (idx_n>=0 && idx_n<(int)f.size()) row["N"]=f[idx_n];

        fout.write_line(FE.format_line(spec, row));
    }
}