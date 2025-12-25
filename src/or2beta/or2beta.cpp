#include "or2beta/or2beta.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/gwasQC.hpp"
#include "utils/FormatEngine.hpp"
#include "utils/StatFunc.hpp"

#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <string_view>
#include <vector>
#include <cstdlib>   // strtod
#include <cerrno>    // errno
#include <limits>

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


void run_or2beta(const Args_Or2Beta& P){
    LineReader lr(P.gwas_file);
    string line;

    if (!lr.getline(line)) {
        LOG_ERROR("Empty GWAS summary file in or2beta.");
        exit(1);
    }
    strip_cr_inplace(line);   
    auto header = split(line);

    // header check
    int idx_snp = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for or2beta.");

    int idx_A1  = find_col(header, P.g_A1);
    require(idx_A1 >= 0,  "GWAS missing required column [" + P.g_A1 + "] for or2beta.");

    int idx_A2  = find_col(header, P.g_A2);
    require(idx_A2 >= 0,  "GWAS missing required column [" + P.g_A2 + "] for or2beta.");

    int idx_or  = find_col(header, P.col_or);
    require(idx_or >= 0, "GWAS missing required column [" + P.col_or + "] for or2beta.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for or2beta.");

    int idx_se   = find_col(header, P.col_se);
    int idx_p    = find_col(header, P.g_p);
    require(idx_se >= 0 || idx_p >= 0,
            "or2beta requires either SE column [" + P.col_se +
            "] or P column [" + P.g_p + "].");

    int idx_n    = find_col(header, P.col_n);

    // ---------- read lines ----------
    std::vector<std::string> lines;
    lines.reserve(1 << 20);

    while (lr.getline(line)){
        if (line.empty()) continue;
        strip_cr_inplace(line); 
        lines.push_back(line);
    }

    size_t n = lines.size();
    LOG_INFO("Loaded " + to_string(n) + " GWAS lines for or2beta.");

    // ======================= QC =======================
    std::vector<bool> keep(n, true);

    bool can_qc = (idx_freq>=0 || idx_p>=0 || idx_n>=0);
    if (can_qc) {
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(
            lines,
            header,
            -1,        // 不 QC beta
            idx_se,    // se 可选
            idx_freq,  // freq 可选
            idx_p,     // p 可选
            idx_n,     // N 可选
            keep,
            P.maf_threshold
        );
    } else {
        LOG_WARN("QC not fully applied: missing freq/p/N columns.");
    }

    // remove dup SNP：用 SNP 作为 key
    if (P.remove_dup_snp) {
        vector<string> snp_vec(n);

        // 不 split，只扫 SNP 列
        int stop = idx_snp;
        std::vector<int> col2slot(stop + 1, -1);
        col2slot[idx_snp] = 0;

        for (size_t i=0;i<n;i++){
            if (!keep[i]) continue;
            std::string_view outs[1] = {};
            int cols = scan_to_stop_col(std::string_view(lines[i]), stop, col2slot, outs, 1);
            if (cols < stop + 1) continue;
            auto v = trim_ws(outs[0]);
            if (!v.empty()) snp_vec[i].assign(v.data(), v.size());
        }
        gwas_remove_dup(lines, header, idx_p, snp_vec, keep);
    }

    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);
    Writer fout(P.out_file, P.format);

    if (!fout.good()) {
        LOG_ERROR("Cannot open output file: " + P.out_file);
        exit(1);
    }

    // writer header
    if (P.format == "gwas") {
        // raw header
        string h;
        for (size_t i=0; i<header.size(); i++){
            if (i) h += "\t";
            h += header[i];
        }
        fout.write_line(h);
    } else {
        string h;
        for (size_t i=0; i<spec.cols.size(); i++){
            if (i) h += "\t";
            h += spec.cols[i];
        }
        fout.write_line(h);
    }

    // 计算需要扫描到的最大列（避免 split）
    int stop = idx_snp;
    stop = std::max(stop, idx_A1);
    stop = std::max(stop, idx_A2);
    stop = std::max(stop, idx_or);
    stop = std::max(stop, idx_freq);
    if (idx_se >= 0) stop = std::max(stop, idx_se);
    if (idx_p  >= 0) stop = std::max(stop, idx_p);
    if (idx_n  >= 0) stop = std::max(stop, idx_n);

    // slot: 0=SNP,1=A1,2=A2,3=OR,4=FREQ,5=SE,6=P,7=N
    std::vector<int> col2slot(stop + 1, -1);
    col2slot[idx_snp]  = 0;
    col2slot[idx_A1]   = 1;
    col2slot[idx_A2]   = 2;
    col2slot[idx_or]   = 3;
    col2slot[idx_freq] = 4;
    if (idx_se >= 0) col2slot[idx_se] = 5;
    if (idx_p  >= 0) col2slot[idx_p]  = 6;
    if (idx_n  >= 0) col2slot[idx_n]  = 7;


    // process lines
    for (size_t i=0; i<n; i++) {
        if (!keep[i]) continue;

        const std::string &ln = lines[i];

        std::string_view outs[8] = {};
        int cols = scan_to_stop_col(std::string_view(ln), stop, col2slot, outs, 8);

        //不再要求 f.size()==header.size()；只要关键列存在即可，避免不必要丢行
        if (cols < stop + 1) continue;

        auto vSNP = trim_ws(outs[0]);
        if (vSNP.empty()) continue;

        if (P.format == "gwas"){
            // 与原逻辑一致：gwas 模式不改行内容（仅过滤不合法行）
            fout.write_line(ln);             // 不 split，直接写
            continue;
        }

        // OR -> beta
        double ORv = NAN;
        if (!parse_double_strict(outs[3], ORv)) continue;   // 用 strtod strict 替代 stod
        if (!(ORv > 0.0) || !std::isfinite(ORv)) continue;

        double beta = std::log(ORv);
        // 计算 se
        double se = NAN;
        if (idx_se >= 0){
            double sev = NAN;
            if (parse_double_strict(outs[5], sev) && sev > 0.0) {
                se = sev;
            }
        }

        if (!std::isfinite(se)) {
            if (idx_p >= 0) {
                double pval = NAN;
                if (parse_double_strict(outs[6], pval) && pval > 0.0 && pval <= 1.0) {
                    double z = StatFunc::p2z_two_tailed(pval);
                    se = (z > 0 ? std::fabs(beta) / z : 999.0);
                } else {
                    se = 999.0;
                }
            } else {
                se = 999.0;
            }
        }

        std::string beta_str = std::to_string(beta);
        std::string se_str   = std::to_string(se);

        FormatEngine::RowView row;                     // FormatEngine
        row.SNP  = {vSNP, true};
        row.A1   = {trim_ws(outs[1]), true};
        row.A2   = {trim_ws(outs[2]), true};

        row.freq = {trim_ws(outs[4]), true};

        // 可选列：N / p
        if (idx_n >= 0) row.N = {trim_ws(outs[7]), true};
        else           row.N = {{}, false};

        if (idx_p >= 0) row.p = {trim_ws(outs[6]), true};
        else           row.p = {{}, false};

        row.beta = {std::string_view(beta_str), true};
        row.se   = {std::string_view(se_str),   true};

        fout.write_line(FE.format_line_fast(spec, row)); // fast path
    }
}