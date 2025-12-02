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

using namespace std;

void process_rsidImpu(const Args_RsidImpu& P,
                    const DBMap& mapdb)
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
    
    vector<bool>   keep(n,false);       // 是否最终进入主输出
    vector<string> rsid_vec(n);         // 匹配到的 rsID

    //================  匹配 dbSNP，确定哪些行有 rsid =================
    for (size_t i=0; i<n; i++){
        if (!keep_qc[i]) {
            keep[i] = false;
            continue;
        }
        
        auto f = split(gwas_lines[i]);

        int max_col = gCHR;
        if (gPOS > max_col) max_col = gPOS;
        if (gA1  > max_col) max_col = gA1;
        if (gA2  > max_col) max_col = gA2;

        if ((int)f.size() <= max_col) {
            keep[i] = false;
            continue;
        }

        // chr-bucket DBMap
        // string key = make_key(f[gCHR], f[gPOS], f[gA1], f[gA2]);
        string chr = norm_chr(f[gCHR]);
        string pos = f[gPOS];
        auto canon = canonical_alleles(f[gA1], f[gA2]);
        string key = pos + ":" + canon.first + ":" + canon.second;

        auto it_chr = mapdb.find(chr);
        if (it_chr != mapdb.end()) {
            auto it = it_chr->second.find(key);
            if (it != it_chr->second.end()) {
                keep[i] = true;             // 只有匹配到 rsID 的才保留
                rsid_vec[i] = it->second;
            } else {
                keep[i] = false;
            }
        } else {
            keep[i] = false;
        }
    }
    
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
        h += "\tSNP";
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
            fout.write_line(gwas_lines[i] + "\t" + rsid_vec[i]);
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