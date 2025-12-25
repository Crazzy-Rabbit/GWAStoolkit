#include "convert/convert.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/gwasQC.hpp"
#include "utils/FormatEngine.hpp"

#include <unordered_map>
#include <string>
#include <deque>
#include <algorithm>
#include <string_view>
#include <vector>
#include <cstdlib>   // strtod
#include <cerrno>    // errno
#include <cmath>     // isfinite, log, fabs
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


void run_convert(const Args_Convert& P){
    LineReader lr(P.gwas_file);
    string line;

    // read header
    if (!lr.getline(line)){
        LOG_ERROR("Empty GWAS summary file in convert.");
        exit(1);
    }
    strip_cr_inplace(line);     //读入时就去 '\r'
    auto header = split(line);
    
    // header check
    int idx_snp  = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for convert.");

    int idx_A1   = find_col(header, P.g_A1);
    require(idx_A1 >= 0, "GWAS missing required column [" + P.g_A1 + "] for convert.");

    int idx_A2   = find_col(header, P.g_A2);
    require(idx_A2 >= 0, "GWAS missing required column [" + P.g_A2 + "] for convert.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for convert.");

    int idx_beta = find_col(header, P.col_beta);
    require(idx_beta >= 0, "GWAS missing required column [" + P.col_beta + "] for convert.");

    int idx_se   = find_col(header, P.col_se);
    require(idx_se >= 0, "GWAS missing required column [" + P.col_se + "] for convert.");

    int idx_p    = find_col(header, P.g_p);
    require(idx_p >= 0, "GWAS missing required column [" + P.g_p + "] for convert.");

    int idx_n    = find_col(header, P.col_n);
    require(idx_n >= 0, "GWAS missing required column [" + P.col_n + "] for convert.");

    // read all line
    std::vector<std::string> lines;
    line.reserve(1 << 20);           // 减少扩容
    while (lr.getline(line)) {
        if (line.empty()) continue;
        strip_cr_inplace(line);    // 读入时去 '\r'
        lines.push_back(line);
    }
    size_t n = lines.size();
    LOG_INFO("Loaded GWAS lines for convert: " + to_string(n));

    // QC, default maf = 0.01
    std::vector<bool> keep(n, true);

    // 如果列都存在，就执行QC，否则给warning
    bool can_qc = (idx_beta>=0 || idx_se>=0 || idx_freq>=0 || idx_p>=0 || idx_n>=0);
    if (can_qc) {
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(lines, header,
                    idx_beta, idx_se, idx_freq, idx_p, idx_n,
                    keep, P.maf_threshold);
    } else {
        LOG_WARN("Cannot perform full QC in convert (missing beta/se/freq/N/P columns).");
    }

    if (P.remove_dup_snp) {
        std::vector<std::string> snp_vec(n);

        // 构建 snp_vec 不再 split，每行只扫到 SNP 列
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

    // out format
    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);
    Writer fout(P.out_file, P.format);

    if (!fout.good()){
        LOG_ERROR("Cannot open output file: " + P.out_file);
        exit(1);
    }

    // writer header
    if (P.format == "gwas") {
        // raw header
        std::string h;
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

    // write context
    if (P.format == "gwas") {
        // gwas 格式直接写原行（不 split）
        for (size_t i=0; i<n; i++){
            if (!keep[i]) continue;
            fout.write_line(lines[i]);
        }
        LOG_INFO("convert finished (format=" + P.format + ").");
        return;
    }

    // 非 gwas：单次扫描取必需列 + FormatEngine fast path
    int stop = std::max({idx_snp, idx_A1, idx_A2, idx_freq, idx_beta, idx_se, idx_p, idx_n});

    std::vector<int> col2slot(stop + 1, -1);
    col2slot[idx_snp]  = 0;
    col2slot[idx_A1]   = 1;
    col2slot[idx_A2]   = 2;
    col2slot[idx_freq] = 3;
    col2slot[idx_beta] = 4;
    col2slot[idx_se]   = 5;
    col2slot[idx_p]    = 6;
    col2slot[idx_n]    = 7;

    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;

        std::string_view outs[8] = {};
        int cols = scan_to_stop_col(std::string_view(lines[i]), stop, col2slot, outs, 8);

        // 不再用 f.size()==header.size()（会导致多列/少列全丢）
        // 只要“至少有我们需要的列”即可；行截断则跳过，避免错位风险。
        if (cols < stop + 1) continue;

        auto vSNP = trim_ws(outs[0]);
        if (vSNP.empty()) continue;

        FormatEngine::RowView row;                      // 新版 FormatEngine
        row.SNP  = {vSNP, true};
        row.A1   = {trim_ws(outs[1]), true};
        row.A2   = {trim_ws(outs[2]), true};
        row.freq = {trim_ws(outs[3]), true};
        row.beta = {trim_ws(outs[4]), true};
        row.se   = {trim_ws(outs[5]), true};
        row.p    = {trim_ws(outs[6]), true};
        row.N    = {trim_ws(outs[7]), true};

        fout.write_line(FE.format_line_fast(spec, row)); // fast path
    }

    LOG_INFO("convert finished (format=" + P.format + ").");
}