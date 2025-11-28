//
//  gwas.cpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/linereader.hpp"
#include "utils/writer.hpp"  // support .gz format out
#include "utils/log.hpp"
#include "utils/util.hpp"
#include "utils/gwasQC.hpp" // basic QC
#include "utils/FormatEngine.hpp"
#include "rsidImpu/gwas.hpp"
#include "rsidImpu/allele.hpp"
#include "rsidImpu/dbsnp.hpp"


#include <fstream>
#include <iostream>
#include <algorithm>
#include <deque>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void process_gwas(const Args_RsidImpu& P,
                    const DBMap& mapdb)
{
    // threads
    #ifdef _OPENMP
    if (P.threads > 0) {
        omp_set_num_threads(P.threads);
        LOG_INFO("Using threads = " + std::to_string(P.threads));
    }
    #endif

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
    int gCHR = find_col(header, P.g_chr);
    int gPOS = find_col(header, P.g_pos);
    int gA1  = find_col(header, P.g_A1);
    int gA2  = find_col(header, P.g_A2);
    int gP   = find_col(header, P.g_p);

    if (gCHR<0||gPOS<0||gA1<0||gA2<0||gP<0){
        LOG_ERROR("GWAS header incomplete.");
        exit(1);
    }

    // format = cojo
    int colFreq = -1, colBeta = -1, colSe = -1, colN = -1;
    bool is_cojo = (P.format == "cojo");

    if (is_cojo){
        colFreq = find_col(header, P.col_freq);
        colBeta = find_col(header, P.col_beta);
        colSe   = find_col(header, P.col_se);
        colN    = find_col(header, P.col_n);

        if (colFreq<0 || colBeta<0 || colSe<0 || colN<0) {
            LOG_ERROR("Error: COJO format requires freq, beta, se, n columns.");
            exit(1);
        }
    }

    // format = smr
    // int colFreq = -1, colBeta = -1, colSe = -1, colN = -1;
    // bool is_smr = (P.format == "smr");

    // if (is_smr){
    //     colFreq = find_col(header, P.col_freq);
    //     colBeta = find_col(header, P.col_beta);
    //     colSe   = find_col(header, P.col_se);
    //     colN    = find_col(header, P.col_n);

    //     if (colFreq<0 || colBeta<0 || colSe<0 || colN<0) {
            // LOG_ERROR("Error: SMR format requires freq, beta, se, n columns.");
            // exit(1);
    //     }
    // }
    
    //================ 2. 读入 GWAS 数据行 =================
    while (reader.getline(line)){
        if (line.empty()) continue;
        gwas_lines.push_back(line);
    }
    
    size_t n = gwas_lines.size();
    LOG_INFO("Loaded GWAS lines (data): " + std::to_string(n));

    //================ 基础 QC：过滤无效 N/beta/se/freq/P =================
    vector<bool> keep_qc(gwas_lines.size(), true);

    int idx_beta = find_col(header, P.col_beta);
    int idx_se   = find_col(header, P.col_se);
    int idx_freq = find_col(header, P.col_freq);
    int idx_pv   = find_col(header, P.g_p);     // pval
    int idx_n    = find_col(header, P.col_n);
    double maf   = P.maf_threshold;
    
    gwas_basic_qc(gwas_lines, header, idx_beta, idx_se, idx_freq, idx_pv, idx_n, keep_qc, maf);

    vector<bool>    keep(n,false);
    vector<string>  rsid_vec(n);
    
    // ---------- 把 QC 过滤结果同步到 keep ----------
    for (size_t i = 0; i < n; i++){
        if (!keep_qc[i]) keep[i] = false;
    }

    //================ 3. 匹配 dbSNP，确定哪些行有 rsid =================
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (long i=0; i<(long)n; i++){
        auto f = split(gwas_lines[i]);
        
        int max_col = gCHR;
        if (gPOS > max_col) max_col = gPOS;
        if (gA1  > max_col) max_col = gA1;
        if (gA2  > max_col) max_col = gA2;
        if (gP   > max_col) max_col = gP;

        if ((int)f.size() <= max_col) continue;

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
                keep[i] = true;
                rsid_vec[i] = it->second;
            }
        }
    }
    
    if (P.remove_dup_snp) {
        gwas_remove_dup(
            gwas_lines,
            header,
            gP,          // P 列 index
            rsid_vec,
            keep
        );
    }
    
    //================ 4. 构建 Writer（自动 txt / gz） =================
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
        string header_line;
        for (size_t j=0; j<header.size(); j++) {
            if (j) header_line += "\t";
            header_line += header[j];
        }
        header_line += "\tSNP";
        fout.write_line(header_line);
    } else {
        // 格式化 header
        string header_line;
        for (size_t j=0; j<spec.cols.size(); j++) {
            if (j) header_line += "\t";
            header_line += spec.cols[j];
        }
        fout.write_line(header_line);
    }

    // outfile 
    for (size_t i=0; i<n; i++){
        gwas_lines[i].erase(
            std::remove(gwas_lines[i].begin(), gwas_lines[i].end(), '\r'),
            gwas_lines[i].end()
        );

        if (!keep[i]) {
            funm.write_line(gwas_lines[i]);
            continue;
        }
        gwas_lines[i].erase(
            std::remove(gwas_lines[i].begin(), gwas_lines[i].end(), '\r'),
            gwas_lines[i].end()
        );
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
        row["freq"] = f[colFreq];
        row["beta"] = f[colBeta];
        row["se"]   = f[colSe];
        row["p"]    = f[gP];
        row["P"]    = f[gP];
        row["N"]    = f[colN];   

        fout.write_line(FE.format_line(spec, row));
    }
}