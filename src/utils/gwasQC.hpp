//
//  qc.hpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_GWASQC_HPP
#define RSIDIMPU_GWASQC_HPP

#include "utils/util.hpp"

#include <string>
#include <deque>
#include <vector>

// =======================================================
// basic QC for gwas summary statistic
// - remove rows where: N, beta, se, freq not finite
// - p must be in [0,1]
// - freq must pass maf threshold if provided
// =======================================================

// ---------------------------
// [API] deque<string> versions (backward compatible)
// ---------------------------
void gwas_basic_qc(
    std::deque<std::string>& lines,
    const std::vector<std::string>& header,
    int idx_beta,
    int idx_se,
    int idx_freq,
    int idx_p,
    int idx_n,
    std::vector<bool>& keep,
    double maf_threshold
);

void gwas_remove_dup(
    std::deque<std::string> &lines,
    const std::vector<std::string> &header,
    int idx_p,                            // gP, can be -1
    std::vector<std::string> &rsid_vec,
    std::vector<bool> &keep
);

// ---------------------------
// [MOD-1] NEW: vector<string> overloads (for faster pipelines)
// 前面已经把多处容器换成 vector<string>，
// 这两个 overload 避免反复 deque<->vector 转换。
// ---------------------------
void gwas_basic_qc(
    std::vector<std::string>& lines,
    const std::vector<std::string>& header,
    int idx_beta,
    int idx_se,
    int idx_freq,
    int idx_p,
    int idx_n,
    std::vector<bool>& keep,
    double maf_threshold
);

void gwas_remove_dup(
    std::vector<std::string> &lines,
    const std::vector<std::string> &header,
    int idx_p,                            // gP, can be -1
    std::vector<std::string> &rsid_vec,
    std::vector<bool> &keep
);

#endif