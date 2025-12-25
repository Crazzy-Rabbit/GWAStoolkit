#include "computeNeff.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/FormatEngine.hpp"
#include "utils/gadgets.hpp"
#include "utils/gwasQC.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <deque>
#include <string_view>
#include <limits>
#include <cstdlib>   // strtod
#include <cerrno>    // errno

using namespace std;
// =======================================================
// [OPT-SHARED] Fast tab scan helpers (no split, no alloc)
// =======================================================

static inline void strip_cr_inplace(std::string &s){
    // [OPT-SHARED-1] 常见情况是行尾 '\r'，O(1) 处理
    if (!s.empty() && s.back() == '\r') { s.pop_back(); return; }
    s.erase(std::remove(s.begin(), s.end(), '\r'), s.end());
}

static inline std::string_view trim_ws(std::string_view sv){
    while (!sv.empty() && (sv.front()==' ' || sv.front()=='\t')) sv.remove_prefix(1);
    while (!sv.empty() && (sv.back() ==' ' || sv.back() =='\t' || sv.back()=='\r')) sv.remove_suffix(1);
    return sv;
}

// 扫描到 stop_col（包含 stop_col），用 col2slot 映射把目标列写入 outs[slot]。
// 返回：扫描到的列数（若行列数>=stop_col+1，则返回 stop_col+1；否则返回实际列数）。
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

            // 够用就停（避免扫描整行）
            if (col - 1 == stop_col) break;
        }
    }
    return col;
}

// [OPT-SHARED-2] 严格 double 解析（不抛异常，比 stod 稳；并做“整串消费”校验）
static inline bool parse_double_strict(std::string_view sv, double &out){
    sv = trim_ws(sv);
    if (sv.empty()) return false;

    errno = 0;
    char *end = nullptr;
    out = std::strtod(sv.data(), &end);

    if (end == sv.data()) return false; // 无法解析
    if ((const char*)end != sv.data() + sv.size()) return false; // [FIX] 必须整串消费，避免 "1abc"
    if (errno == ERANGE) return false;
    if (!std::isfinite(out)) return false;
    return true;
}

// [OPT-SHARED-3] 预计算某列的 [start,len]（用于原地替换，不 split）
// 成功返回 true；失败返回 false（行列不够）
static inline bool get_col_span(std::string_view line, int col_idx, uint32_t &st, uint32_t &len){
    if (col_idx < 0) return false;
    size_t start = 0;
    for (int c=0; c<col_idx; ++c){
        size_t p = line.find('\t', start);
        if (p == std::string_view::npos) return false;
        start = p + 1;
    }
    size_t end = line.find('\t', start);
    if (end == std::string_view::npos) end = line.size();
    st  = (uint32_t)start;
    len = (uint32_t)(end - start);
    return true;
}


// Neff
static inline double calc_neff(double cs, double ct){
    double s = cs + ct;
    if (s <= 0) return 0;
    return 4.0 * cs * ct / s;
}

// standardization of beta and se
static inline bool std_effect(
    double freq, double beta_old, double se_old, double Neff,
    double &beta_new, double &se_new
) {
    if (freq <= 0.0 || freq >= 1.0) return false;
    if (se_old <= 0.0) return false;
    if (!std::isfinite(Neff) || Neff <= 0.0) return false;

    double z = beta_old / se_old;
    double denom = 2.0 * freq * (1.0-freq) * (Neff + z*z);
    if (denom <= 0.0) return false;

    se_new   = 1.0 / sqrt(denom);
    beta_new = z * se_new;
    return true;
}

void run_computeNeff(const Args_CalNeff& P)
{

    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);
    Writer fout(P.out_file, P.format);

    if (!fout.good()){
        LOG_ERROR("Cannot open output: " + P.out_file);
        return;
    }

    // header
    LineReader reader(P.gwas_file);
    std::string line;
    if (!reader.getline(line)) {
        LOG_ERROR("Empty GWAS file: " + P.gwas_file);
        exit(1);
    }
    strip_cr_inplace(line);
    auto header = split(line);

    bool has_N = false;
    int idx_N  = -1;

    for (size_t i = 0; i < header.size(); i++){
        std::string hlow = header[i];
        std::transform(hlow.begin(), hlow.end(), hlow.begin(), ::tolower);

        if (hlow == "n"){
            has_N = true;
            idx_N = i;
            break;
        }
    }

    // header check
    int idx_snp = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for computeNeff.");

    int idx_A1 = find_col(header, P.g_A1);
    require(idx_A1 >= 0, "GWAS missing required column [" + P.g_A1 + "] for computeNeff.");

    int idx_A2 = find_col(header, P.g_A2);
    require(idx_A2 >= 0, "GWAS missing required column [" + P.g_A2 + "] for computeNeff.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for computeNeff.");

    int idx_beta = find_col(header, P.col_beta);
    require(idx_beta >= 0, "GWAS missing required column [" + P.col_beta + "] for computeNeff.");

    int idx_se = find_col(header, P.col_se);
    require(idx_se >= 0, "GWAS missing required column [" + P.col_se + "] for computeNeff.");

    int idx_p = find_col(header, P.g_p);
    require(idx_p >= 0, "GWAS missing required column [" + P.g_p + "] for computeNeff.");

    int idx_case    = -1;
    int idx_control = -1;

    // 读入所有数据行
    std::vector<std::string> lines;
    lines.reserve(1 << 20);
    std::vector<std::pair<uint32_t,uint32_t>> n_span;

    if (P.format == "gwas" && has_N) n_span.reserve(1 << 20);

    while (reader.getline(line)){
        if (line.empty()) continue;
        strip_cr_inplace(line);
        lines.push_back(line);

        // 预计算 N 列 span（避免后面 split 重建整行）
        if (P.format == "gwas" && has_N) {
            uint32_t st=0, len=0;
            bool ok = get_col_span(std::string_view(lines.back()), idx_N, st, len);
            if (!ok) { st = std::numeric_limits<uint32_t>::max(); len = 0; }
            n_span.emplace_back(st, len);
        }
    }
    size_t n = lines.size();
    LOG_INFO("Loaded " + to_string(n) + " GWAS lines for computeNeff.");

    // -------------------------
    // Case 1: per-SNP NEFF 计算
    // -------------------------
    if (P.is_column){
        idx_case    = find_col(header, P.case_col);
        idx_control = find_col(header, P.control_col);
        if (idx_case < 0) {
            LOG_ERROR("Cannot find case column: " + P.case_col);
            exit(1);
        }
        if (idx_control < 0) {
            LOG_ERROR("Cannot find control column: " + P.control_col);
            exit(1);
        }
    }

    // -------------------------
    // Case 2: fixed case/control
    // -------------------------
    double Neff_fixed = std::numeric_limits<double>::quiet_NaN();
    if (P.is_single){
        Neff_fixed = calc_neff(P.case_n, P.control_n);
        if (!std::isfinite(Neff_fixed) || Neff_fixed <= 0.0){
            LOG_ERROR("Invalid fixed case/control: " + 
                std::to_string(P.case_n) + "," + std::to_string(P.control_n));
            exit(1);
        }
        LOG_INFO("Fixed-mode Neff = " + std::to_string(Neff_fixed));

    }

    //================ QC + 去重 (按 SNP) ================
    std::vector<bool> keep(n, true);
    
    bool can_qc = (idx_beta>=0 || idx_se>=0 || idx_freq>=0 || idx_p>=0);
    if (can_qc){
        int idx_n_safe = (has_N ? idx_N : -1); 
        
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(lines, header,
                    idx_beta, idx_se, idx_freq, idx_p, idx_n_safe,
                    keep, P.maf_threshold);
    } else {
        LOG_WARN("Cannot perform full QC in computeNeff (missing beta/se/freq/N/p columns).");
    }

    // remove dup SNP（用 SNP 作为 key）
    if (P.remove_dup_snp){
        std::vector<string> snp_vec(n);

        // snp_vec 构建：不 split，只扫 SNP 列
        int stop = idx_snp;
        std::vector<int> col2slot(stop + 1, -1);
        col2slot[idx_snp] = 0;

        for (size_t i=0; i<n; i++){
            if (!keep[i]) continue;
            std::string_view outs[1] = {};
            int cols = scan_to_stop_col(std::string_view(lines[i]), stop, col2slot, outs, 1);
            if (cols < stop + 1) continue;
            auto v = trim_ws(outs[0]);
            if (!v.empty()) snp_vec[i].assign(v.data(), v.size());
        }
        gwas_remove_dup(lines, header, idx_p, snp_vec, keep);
    }


    // writer header
    if (P.format == "gwas") {
        // raw header
        string h;
        for (size_t i=0; i<header.size(); i++){
            if (i) h += "\t";
            h += header[i];
        }

        if (!has_N) h += "\tN";
        fout.write_line(h);
    } else {
        string h;
        for (size_t i=0; i<spec.cols.size(); i++){
            if (i) h += "\t";
            h += spec.cols[i];
        }
        fout.write_line(h);
    }

    // process GWAS per row
    // 分支：gwas 输出只需要 SNP + case/control（无需扫 A1/A2/freq/beta/se/p）
    int stop_min = idx_snp;
    if (P.is_column) stop_min = std::max({stop_min, idx_case, idx_control});

    int stop_full = std::max({idx_snp, idx_A1, idx_A2, idx_freq, idx_beta, idx_se, idx_p});
    if (P.is_column) stop_full = std::max({stop_full, idx_case, idx_control});

    // 为两种 stop 准备 col2slot（减少每行构造开销）
    // --- minimal slots: SNP, CASE, CONTROL ---
    std::vector<int> col2slot_min(stop_min + 1, -1);
    // slot: 0=SNP, 1=CASE, 2=CONTROL
    col2slot_min[idx_snp] = 0;
    if (P.is_column){
        col2slot_min[idx_case] = 1;
        col2slot_min[idx_control] = 2;
    }

    // --- full slots: SNP,A1,A2,freq,beta,se,p,CASE,CONTROL ---
    std::vector<int> col2slot_full(stop_full + 1, -1);
    // slot: 0=SNP,1=A1,2=A2,3=freq,4=beta,5=se,6=p,7=CASE,8=CONTROL
    col2slot_full[idx_snp]  = 0;
    col2slot_full[idx_A1]   = 1;
    col2slot_full[idx_A2]   = 2;
    col2slot_full[idx_freq] = 3;
    col2slot_full[idx_beta] = 4;
    col2slot_full[idx_se]   = 5;
    col2slot_full[idx_p]    = 6;
    if (P.is_column){
        col2slot_full[idx_case] = 7;
        col2slot_full[idx_control] = 8;
    }

    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;

        const std::string &ln = lines[i];
        
        // 计算当前 SNP 的 Neff
        double Neff = NAN;
        if (P.is_single){
            Neff = Neff_fixed;
        } else if (P.is_column){
            std::string_view outs_min[3] = {};
            int cols = scan_to_stop_col(std::string_view(ln), stop_min, col2slot_min, outs_min, 3);
            if (cols < stop_min + 1) continue;

            double cs=0.0, ct=0.0;
            if (!parse_double_strict(outs_min[1], cs)) continue;
            if (!parse_double_strict(outs_min[2], ct)) continue;
            Neff = calc_neff(cs, ct);
        }

        if (!std::isfinite(Neff) || Neff <= 0.0) continue;
        
        // ----------- gwas 输出：原地替换/追加 N（不 split） -----------
        if (P.format == "gwas"){
            std::string neff_str = std::to_string(Neff);

            if (has_N){
                // 用预计算 span 直接替换 N 列
                if (i < n_span.size()){
                    auto [st, len] = n_span[i];
                    if (st != std::numeric_limits<uint32_t>::max()){
                        std::string out = ln; // 保持 lines 不被破坏（可改为就地改写）
                        out.replace((size_t)st, (size_t)len, neff_str);
                        fout.write_line(out);
                    } else {
                        // 行截断导致找不到 N 列：安全兜底 -> 追加
                        fout.write_line(ln + "\t" + neff_str); // [FIX-NEFF-2]
                    }
                } else {
                    fout.write_line(ln + "\t" + neff_str);
                }
            } else {
                // 原来没有 N -> 追加
                fout.write_line(ln + "\t" + neff_str);
            }
            continue;
        }

        // ----------- 非 gwas 输出：需要更多字段（SNP/A1/A2/freq/beta/se/p） -----------
        std::string_view outs[9] = {};
        int cols = scan_to_stop_col(std::string_view(ln), stop_full, col2slot_full, outs, 9);
        if (cols < stop_full + 1) continue;

        auto vSNP = trim_ws(outs[0]);
        if (vSNP.empty()) continue;

        // 标准化 beta/se（若失败则保留旧值）
        double freq_old=0, beta_old=0, se_old=0;
        if (!parse_double_strict(outs[3], freq_old)) continue;
        if (!parse_double_strict(outs[4], beta_old)) continue;
        if (!parse_double_strict(outs[5], se_old))   continue;

        std::string beta_str, se_str, neff_str;
        neff_str = std::to_string(Neff);

        double beta_new=0, se_new=0;
        bool ok_std = std_effect(freq_old, beta_old, se_old, Neff, beta_new, se_new);

        if (ok_std){
            beta_str = std::to_string(beta_new);
            se_str   = std::to_string(se_new);
        }

        FormatEngine::RowView row;                       // 新版 FormatEngine
        row.SNP  = {vSNP, true};
        row.A1   = {trim_ws(outs[1]), true};
        row.A2   = {trim_ws(outs[2]), true};
        row.freq = {trim_ws(outs[3]), true};

        if (ok_std){
            row.beta = {std::string_view(beta_str), true};
            row.se   = {std::string_view(se_str),   true};
        } else {
            row.beta = {trim_ws(outs[4]), true};
            row.se   = {trim_ws(outs[5]), true};
        }

        row.p = {trim_ws(outs[6]), true};
        row.N = {std::string_view(neff_str), true};

        fout.write_line(FE.format_line_fast(spec, row));  // fast path
    }
}