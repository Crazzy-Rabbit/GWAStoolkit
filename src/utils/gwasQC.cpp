//
//  gwasQC.cpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/gwasQC.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"

#include <algorithm>
#include <cmath>
#include <cerrno>      // [MOD-INC] errno / ERANGE
#include <cstdlib>     // [MOD-INC] strtod
#include <cstring>     // [MOD-INC] memcpy
#include <string_view> // [MOD-INC] string_view
#include <unordered_map>
#include <vector>      // [MOD-INC] 支持 vector<string> overload
#include <utility>

using namespace std;

// =======================================================
// [OPT-SHARED] Fast parsing helpers (no split, no exceptions)
// =======================================================

static inline std::string_view trim_ws(std::string_view sv){
    while (!sv.empty() && (sv.front()==' ' || sv.front()=='\t')) sv.remove_prefix(1);
    while (!sv.empty() && (sv.back() ==' ' || sv.back() =='\t' || sv.back()=='\r')) sv.remove_suffix(1);
    return sv;
}

// strtod 严格解析 double：不抛异常、速度远快于 stod+try/catch
static inline bool parse_double_strict(std::string_view sv, double &out){
    sv = trim_ws(sv);
    if (sv.empty()) return false;

    // 绝大多数数字字段很短：用栈 buffer 避免堆分配
    char buf[128];
    if (sv.size() < sizeof(buf)){
        std::memcpy(buf, sv.data(), sv.size());
        buf[sv.size()] = '\0';

        errno = 0;
        char *end = nullptr;
        out = std::strtod(buf, &end);

        if (end == buf) return false;          // 没读到任何数字
        if (*end != '\0') return false;        // [FIX] 必须整串消费，避免 "1abc"
        if (errno == ERANGE) return false;
        return std::isfinite(out);
    }

    // 极少见：字段特别长才走这里（仍然无异常）
    std::string tmp(sv);
    errno = 0;
    char *end = nullptr;
    out = std::strtod(tmp.c_str(), &end);
    if (end == tmp.c_str()) return false;
    if (*end != '\0') return false;
    if (errno == ERANGE) return false;
    return std::isfinite(out);
}

// 只扫描到 stop_col（包含 stop_col），把目标列写入 outs[slot]。
// col2slot[col] == -1 表示该列不需要。
// 返回：解析到的列数（至少是已扫描到的列数）。
static inline int scan_to_stop_col(
    std::string_view line,
    int stop_col,
    const std::vector<int> &col2slot,
    std::string_view *outs,
    int nouts
){
    (void)nouts; // outs 的长度由调用方保证
    size_t start = 0;
    int col = 0;

    for (size_t j = 0; j <= line.size(); ++j){
        if (j == line.size() || line[j] == '\t'){
            if (col <= stop_col){
                int slot = col2slot[col];
                if (slot >= 0){
                    outs[slot] = std::string_view(line.data() + start, j - start);
                }
            }
            ++col;
            start = j + 1;

            // [OPT] 到 stop_col 就停：避免扫描整行
            if (col - 1 == stop_col) break;
        }
    }
    return col;
}

// =======================================================
// [MOD] Internal templated impl: support deque<string> and vector<string>
// =======================================================

template <class LinesT>
static void gwas_basic_qc_impl(
    LinesT &lines,
    const vector<string> &header,
    int idx_beta,
    int idx_se,
    int idx_freq,
    int idx_p,
    int idx_n,
    vector<bool> &keep,
    double maf_threshold
){
    (void)header; // header 在此函数中不使用，避免 warning

    const size_t n = lines.size();
    size_t kept = 0, dropped = 0;

    // 预计算 stop_col + col2slot，避免每行 split
    int stop = -1;
    if (idx_beta >= 0) stop = std::max(stop, idx_beta);
    if (idx_se   >= 0) stop = std::max(stop, idx_se);
    if (idx_freq >= 0) stop = std::max(stop, idx_freq);
    if (idx_p    >= 0) stop = std::max(stop, idx_p);
    if (idx_n    >= 0) stop = std::max(stop, idx_n);

    // 没有任何可 QC 的列：保持 keep 不变，只统计
    if (stop < 0){
        for (size_t i=0;i<n;++i) if (keep[i]) kept++;
        LOG_INFO("Basic QC done: " + std::to_string(kept) + " passed, 0 removed.");
        return;
    }

    // slot: 0=beta,1=se,2=freq,3=p,4=n
    std::vector<int> col2slot(stop + 1, -1);
    if (idx_beta >= 0) col2slot[idx_beta] = 0;
    if (idx_se   >= 0) col2slot[idx_se]   = 1;
    if (idx_freq >= 0) col2slot[idx_freq] = 2;
    if (idx_p    >= 0) col2slot[idx_p]    = 3;
    if (idx_n    >= 0) col2slot[idx_n]    = 4;


    for (size_t i = 0; i < n; i ++) {
        if (!keep[i]) continue;

        // 不再拷贝 ln，不再 erase/remove('\r')：直接 string_view + trim
        std::string_view lv(lines[i]);

        std::string_view outs[5] = {};
        int cols = scan_to_stop_col(lv, stop, col2slot, outs, 5);

        // [FIX-1] 行列不足（等价于旧 split 后 idx 越界）：QC fail
        if (cols < stop + 1) {
            keep[i] = false;
            dropped++;
            continue;
        }

        bool qc_fail = false;
        double v_beta=0, v_se=0, v_freq=0, v_p=0, v_n=0;

        // [OPT-5] 列不存在(idx<0) → 忽略；列存在 → 严格解析数值
        if (idx_beta >= 0 && !parse_double_strict(outs[0], v_beta)) qc_fail = true;
        if (idx_se   >= 0 && !parse_double_strict(outs[1], v_se))   qc_fail = true;
        if (idx_freq >= 0 && !parse_double_strict(outs[2], v_freq)) qc_fail = true;
        if (idx_p    >= 0 && !parse_double_strict(outs[3], v_p))    qc_fail = true;
        if (idx_n    >= 0 && !parse_double_strict(outs[4], v_n))    qc_fail = true;

        if (qc_fail){
            keep[i] = false;
            dropped++;
            continue;
        }

        // p ∈ [0,1]
        if (idx_p >= 0) {
            if (v_p < 0.0 || v_p > 1.0){
                keep[i] = false; dropped++; continue;
            }
        }

        // MAF: freq ∈ [maf, 1-maf]
        if (idx_freq >= 0) {
            if (v_freq < maf_threshold || v_freq > (1.0 - maf_threshold)) {
                keep[i] = false; dropped++; continue;
            }
        }

        kept++;
    }

    LOG_INFO("Basic QC done: " + std::to_string(kept) + " passed, " + std::to_string(dropped) + " removed.");
}

template <class LinesT>
static void gwas_remove_dup_impl(
    LinesT &lines,
    const vector<string> &header,
    int idx_p,
    vector<string> &rsid_vec,
    vector<bool> &keep
){
    (void)header;

    const size_t n = lines.size();
    size_t dropped = 0;

    // [FIX-2] idx_p < 0：旧实现会访问 f[-1] 崩溃；这里改为“保留首次出现，后续重复删掉”
    if (idx_p < 0){
        struct SVHash { size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); } };
        struct SVEq   { bool   operator()(std::string_view a, std::string_view b) const noexcept { return a==b; } };

        std::unordered_map<std::string_view, size_t, SVHash, SVEq> seen;

        // [OPT-6] reserve 减少 rehash
        size_t active = 0;
        for (size_t i=0;i<n;++i) if (keep[i] && !rsid_vec[i].empty()) ++active;
        seen.reserve(active * 2 + 1);

        for (size_t i=0;i<n;++i){
            if (!keep[i]) continue;
            if (rsid_vec[i].empty()) continue;

            std::string_view snp(rsid_vec[i]);
            auto [it, inserted] = seen.emplace(snp, i);
            if (!inserted){
                keep[i] = false;
                dropped++;
            }
        }

        LOG_INFO("Duplicate SNPs removal done (no P column). Removed = " + std::to_string(dropped));
        return;
    }

    // idx_p >= 0：按最小 p 保留
    struct SVHash { size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); } };
    struct SVEq   { bool   operator()(std::string_view a, std::string_view b) const noexcept { return a==b; } };

    std::unordered_map<std::string_view, std::pair<double, size_t>, SVHash, SVEq> best;

    // [OPT-6] reserve：避免频繁 rehash（大量 SNP 时非常重要）
    size_t active = 0;
    for (size_t i=0;i<n;++i) if (keep[i] && !rsid_vec[i].empty()) ++active;
    best.reserve(active * 2 + 1);

    // [OPT-7] remove_dup 只需要扫描 P 列
    int stop = idx_p;
    std::vector<int> col2slot(stop + 1, -1);
    col2slot[idx_p] = 0;

    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;
        if (rsid_vec[i].empty()) continue;

        std::string_view lv(lines[i]);

        std::string_view outs[1] = {};
        int cols = scan_to_stop_col(lv, stop, col2slot, outs, 1);

        if (cols < stop + 1){
            keep[i] = false;
            dropped++;
            continue;
        }

        double p = 0.0;
        if (!parse_double_strict(outs[0], p)){
            keep[i] = false;
            dropped++;
            continue;
        }

        if (!std::isfinite(p)){
            keep[i] = false;
            dropped++;
            continue;
        }

        std::string_view snp(rsid_vec[i]);

        // [OPT-8] 单次查找/插入：避免 best.count + best[snp] 双重哈希
        auto [it, inserted] = best.emplace(snp, std::make_pair(p, i));
        if (inserted) continue;

        auto &old = it->second;
        if (p < old.first) {
            keep[old.second] = false;
            old = {p, i};
            dropped++;
        } else {
            keep[i] = false;
            dropped++;
        }
    }

    LOG_INFO("Duplicate SNPs removal done. Removed = " + std::to_string(dropped));
}

// =======================================================
// Public API: keep your original signatures + add vector overloads
// =======================================================
// ---------------------------
// basic QC + MAF  (deque<string> version)
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
    // [MOD] 统一走 template 实现
    gwas_basic_qc_impl(lines, header, idx_beta, idx_se, idx_freq, idx_p, idx_n, keep, maf_threshold);
}

// [MOD] 新增：vector<string> 版本（让你前面优化后的模块不必再做 deque<->vector 转换）
// 你需要在 gwasQC.hpp 里也声明这个 overload 才能在别的 cpp 里调用。
void gwas_basic_qc(
    vector<string> &lines,
    const vector<string> &header,
    int idx_beta,
    int idx_se,
    int idx_freq,
    int idx_p,
    int idx_n,
    vector<bool> &keep,
    double maf_threshold
){
    gwas_basic_qc_impl(lines, header, idx_beta, idx_se, idx_freq, idx_p, idx_n, keep, maf_threshold);
}

// ---------------------------
// remove dup SNPs, retain small p  (deque<string> version)
// ---------------------------
void gwas_remove_dup(
    deque<string> &lines,
    const vector<string> &header,
    int idx_p,
    vector<string> &rsid_vec,
    vector<bool> &keep
) {
    // [MOD] 统一走 template 实现
    gwas_remove_dup_impl(lines, header, idx_p, rsid_vec, keep);
}

// [MOD] 新增：vector<string> 版本（同上，避免容器转换）
// 你需要在 gwasQC.hpp 里也声明这个 overload 才能在别的 cpp 里调用。
void gwas_remove_dup(
    vector<string> &lines,
    const vector<string> &header,
    int idx_p,
    vector<string> &rsid_vec,
    vector<bool> &keep
) {
    gwas_remove_dup_impl(lines, header, idx_p, rsid_vec, keep);
}