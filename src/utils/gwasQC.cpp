//
//  gwasQC.cpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/gwasQC.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"

#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cmath>

using namespace std;

static bool parse_num(const string &s, double &v){
    try{
        v = stod(s);
        return std::isfinite(v);
    } catch (...) { return false; }
}

// ---------------------------
// basic QC + MAF 
// ---------------------------
void gwas_basic_qc(
    deque<string> &lines,
    const vector<string> &header,
    int idx_beta,
    int idx_se,
    int idx_freq,
    int idx_p,
    int idx_n,
    vector<bool> &keep,
    double maf_threshold
){
    size_t n = lines.size();
    size_t kept = 0, dropped = 0;

    for (size_t i = 0; i < n; i ++) {
        if (!keep[i]) continue;

        string ln = lines[i];
        ln.erase(remove(ln.begin(), ln.end(), '\r'), ln.end());
        auto f = split(ln);

        int maxCol = idx_beta;
        maxCol = max(maxCol, idx_se);
        maxCol = max(maxCol, idx_freq);
        maxCol = max(maxCol, idx_p);
        maxCol = max(maxCol, idx_n);

        if ((int)f.size() <= maxCol) {
            keep[i] = false; dropped++; continue;
        }

        double v_beta, v_se, v_freq, v_p, v_n;

        bool ok = 
            parse_num(f[idx_beta], v_beta) &&
            parse_num(f[idx_se],   v_se)   &&
            parse_num(f[idx_freq], v_freq) &&
            parse_num(f[idx_p],    v_p)    &&
            parse_num(f[idx_n],    v_n);
        
        // not finite
        if (!ok) { keep[i] = false; dropped++; continue; }
        
        // p (0, 1)
        if (v_p < 0.0 || v_p > 1.0) {
            keep[i] = false; dropped++; continue;
        }

        // MAF filter
        if (v_freq < maf_threshold || v_freq > (1.0 - maf_threshold)) {
            keep[i]=false; dropped++; continue;
        }

        kept++;
    }
    
    LOG_INFO("Basic QC done: " + std::to_string(kept) + " passed, " + std::to_string(dropped) + " removed.");
}

// remove dup SNPs, retain small p 
void gwas_remove_dup(
    deque<string> &lines,
    const vector<string> &header,
    int idx_p,
    vector<string> &rsid_vec,
    vector<bool> &keep
) {
    unordered_map<string, pair<double, size_t>> best;
    size_t n = lines.size();
    size_t dropped = 0;

    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;            // 只看 QC+匹配都通过的行
        if (rsid_vec[i].empty()) continue; // 没有 rsid 的行也跳过

        string ln = lines[i];
        ln.erase(remove(ln.begin(), ln.end(), '\r'), ln.end());
        auto f = split(ln);

        if ((int)f.size() <= idx_p){
            keep[i] = false;
            dropped++;
            continue;
        }

        double p;
        try { p = stod(f[idx_p]); }
        catch(...) { keep[i] = false; dropped++; continue; }

        if (!std::isfinite(p)){
            keep[i] = false;
            dropped++;
            continue;
        }
        
        string snp = rsid_vec[i];

        if (!best.count(snp)){
            best[snp] = {p, i};
        }else{
            auto &old = best[snp];
            if (p < old.first) {
                // 当前行 p 更小 → 保留当前，丢掉之前的
                keep[old.second] = false;
                old = {p, i};
                dropped++;
            } else {
                // 当前行 p 更大 → 丢掉当前
                keep[i] = false;
                dropped++;
            }
        }
    }

    LOG_INFO("Duplicate SNPs removal done. Removed = " + std::to_string(dropped));
}