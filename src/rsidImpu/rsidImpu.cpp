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

#include <algorithm>
#include <cctype>
// #include <charconv>             //  from_chars
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <string>

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

// ---------------- small utilities ----------------
static inline void strip_cr_inplace(std::string &s) {
    // [OPT-0]更快：绝大多数 \r在行尾
    if (!s.empty() && s.back() == '\r') { s.pop_back(); return; }
    s.erase(std::remove(s.begin(), s.end(), '\r'), s.end());
}

static inline vector<string> split_tab(const string& s)
{
    // 仅用于 header（极少次数），保留原实现
    vector<string> out;
    size_t start = 0;
    while (true) {
        size_t pos = s.find('\t', start);
        if (pos == string::npos) {
            out.emplace_back(s.substr(start));
            break;
        }
        out.emplace_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

//  预计算 SNP 列 span：start,len（避免每条输出重复 find '\t'）
static inline bool get_col_span(std::string_view line, int col_idx, uint32_t &st, uint32_t &len){
    if (col_idx < 0) return false;

    size_t start = 0;
    for (int c = 0; c < col_idx; ++c){
        size_t p = line.find('\t', start);
        if (p == std::string_view::npos) return false;
        start = p + 1;
    }

    size_t end = line.find('\t', start);
    if (end == std::string_view::npos) end = line.size();

    // ✅ 正确判断：end 不能小于 start
    if (start > line.size() || end < start) return false;

    st  = static_cast<uint32_t>(start);
    len = static_cast<uint32_t>(end - start);
    return true;
}


static inline char low(char c){
    return (char)std::tolower((unsigned char)c);
}

//  parse int64
static inline std::string_view trim_ws(std::string_view sv){
    while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t')) sv.remove_prefix(1);
    while (!sv.empty() && (sv.back()  == ' ' || sv.back()  == '\t' || sv.back() == '\r')) sv.remove_suffix(1);
    return sv;
}

static inline bool parse_i64(std::string_view sv, int64_t &out) {
    sv = trim_ws(sv);
    if (sv.empty()) return false;

    size_t i = 0;
    bool neg = false;
    if (sv[i] == '+' || sv[i] == '-') {
        neg = (sv[i] == '-');
        ++i;
        if (i == sv.size()) return false;
    }

    int64_t val = 0;
    for (; i < sv.size(); ++i) {
        unsigned char c = (unsigned char)sv[i];
        if (c < '0' || c > '9') return false;
        int digit = int(c - '0');

        // overflow check
        if (val > (LLONG_MAX - digit) / 10) return false;
        val = val * 10 + digit;
    }

    out = neg ? -val : val;
    return true;
}


static inline bool starts_with_ci(std::string_view s, std::string_view p){
    if (s.size() < p.size()) return false;
    for (size_t i=0;i<p.size();++i){
        if (low(s[i]) != low(p[i])) return false;
    }
    return true;
}

static inline bool eq_ci(std::string_view a, std::string_view b){
    if (a.size() != b.size()) return false;
    for (size_t i=0;i<a.size();++i){
        if (low(a[i]) != low(b[i])) return false;
    }
    return true;
}

// 严格整数解析：必须整串都是数字（不允许 "1.11" 这种）
static inline bool parse_int_strict(std::string_view sv, int &out) {
    int64_t v = 0;
    if (!parse_i64(sv, v)) return false;
    if (v < INT_MIN || v > INT_MAX) return false;
    out = (int)v;
    return true;
}

//[OPT-1] 用string_view 快速 canonical_chr （避免 string 分配）
static inline int canonical_chr_code_sv(std::string_view sv) {
    sv = trim_ws(sv);
    if (sv.empty()) return -1;

    // 去掉 CHR 前缀（大小写不敏感）
    if (starts_with_ci(sv, "CHR")) {
        sv.remove_prefix(3);
        sv = trim_ws(sv);
        if (sv.empty()) return -1;
    }

    // RefSeq: NC_000001.11 -> 1; NC_000023.11 -> 23; NC_012920.1 -> 25(MT)
    if (starts_with_ci(sv, "NC_")) {
        sv.remove_prefix(3);
        sv = trim_ws(sv);

        // 至少需要 6 位数字
        if (sv.size() < 6) return -1;

        std::string_view num6 = sv.substr(0, 6);

        int v = 0;
        // ✅ 必须解析 num6，而不是 sv（sv 里会有 ".11"）
        if (!parse_int_strict(num6, v)) return -1;

        // ✅ 注意：这里不能先限制 v<=25，因为 12920 需要映射到 MT
        if (1 <= v && v <= 22) return v;
        if (v == 23) return 23;
        if (v == 24) return 24;
        if (v == 12920) return 25; // NC_012920.* -> MT
        return -1;
    }

    // 常用别名：M / MT / MTDNA
    if (eq_ci(sv, "X")) return 23;
    if (eq_ci(sv, "Y")) return 24;
    if (eq_ci(sv, "M") || eq_ci(sv, "MT") || eq_ci(sv, "MTDNA")) return 25;

    // 数字染色体（严格）
    int v = 0;
    if (!parse_int_strict(sv, v)) return -1;
    if (v <= 0 || v > 25) return -1;
    return v;
}

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
        line.replace(start, line.size() - start, value);
    } else {
        line.replace(start, end - start, value);
    }
}


// =======================================================
// [OPT-1] Tab 扫描器：支持“扫到某一列就停”（dbSNP 核心提速）
// =======================================================
struct TabState {
    size_t i = 0;
    int col  = 0;
};

// 扫描到 stop_col（包含 stop_col），遇到目标列则写入对应 view；支持多次调用继续扫描（靠 TabState）
static inline void scan_upto_col(
    std::string_view line,
    int stop_col,
    int cCHR, int cPOS, int cA1, int cA2, int cRS,
    std::string_view& vCHR,
    std::string_view& vPOS,
    std::string_view& vA1,
    std::string_view& vA2,
    std::string_view& vRS,
    TabState& st
){
    size_t start = st.i;
    for (size_t j = st.i; j <= line.size(); ++j){
        if (j == line.size() || line[j] == '\t'){
            std::string_view sv(line.data() + start, j - start);
            int c = st.col;

            if (c == cCHR) vCHR = sv;
            if (c == cPOS) vPOS = sv;
            if (c == cA1)  vA1  = sv;
            if (c == cA2)  vA2  = sv;
            if (c == cRS)  vRS  = sv;

            st.col++;
            start = j + 1;
            
            if (c == stop_col){
                st.i = start;
                return;
            }
        }
    }
    st.i = line.size();
}

// 额外：输出格式化时需要更多列 -> 扫 7 个目标列
static inline void scan_upto_col7(
    std::string_view line,
    int stop_col,
    int cA1, int cA2, int cFreq, int cBeta, int cSe, int cP, int cN,
    std::string_view& vA1,
    std::string_view& vA2,
    std::string_view& vFreq,
    std::string_view& vBeta,
    std::string_view& vSe,
    std::string_view& vP,
    std::string_view& vN,
    TabState& st
){
    size_t start = st.i;
    for (size_t j = st.i; j <= line.size(); ++j){
        if (j == line.size() || line[j] == '\t') {
            std::string_view sv(line.data() + start, j - start);
            int c = st.col;

            if (c == cA1)   vA1   = sv;
            if (c == cA2)   vA2   = sv;
            if (c == cFreq) vFreq = sv;
            if (c == cBeta) vBeta = sv;
            if (c == cSe)   vSe   = sv;
            if (c == cP)    vP    = sv;
            if (c == cN)    vN    = sv;

            st.col++;
            start = j + 1;

            if (c == stop_col) {
                st.i = start;
                return;
            }
        }
    }
    st.i = line.size();
}

void process_rsidImpu(const Args_RsidImpu& P)
{    
    //================ 1. 读取 GWAS header =================
    LineReader reader(P.gwas_file);
    std::string line;
    if (!reader.getline(line)) {
        LOG_ERROR("Empty GWAS summary file.");
        exit(1);
    }
    strip_cr_inplace(line);

    auto header = split_tab(line);

    // chech exist of SNP col
    bool has_SNP = false;
    int  idx_SNP = -1;

    for (size_t i = 0; i < header.size(); i++){
        string hlow = header[i];
        std::transform(hlow.begin(), hlow.end(), hlow.begin(), ::tolower);

        if (hlow == "snp"){
            has_SNP = true;
            idx_SNP = i;
            break;
        }
    }
    
    // find the col, CHR, POS, A1, A2, P
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
    std::vector<std::string> gwas_lines;
    gwas_lines.reserve(1 << 20); // 可调：减少扩容次数（不影响逻辑）

    // 读 GWAS 时直接构建 gwas_vec，避免第二次 split
    std::vector<GWASRecord> gwas_vec;
    gwas_vec.reserve(1 << 20);

    // 如果需要覆盖 SNP 列，预计算每行 span
    std::vector<std::pair<uint32_t, uint32_t>> snp_span;
    if (P.format == "gwas" && has_SNP) snp_span.reserve(1 << 20);

    while (reader.getline(line)){
        if (line.empty()) continue;
        strip_cr_inplace(line);

        size_t idx = gwas_lines.size();
        gwas_lines.push_back(line);

        // 预存 SNP 列位置
        if (P.format == "gwas" && has_SNP) {
            uint32_t st=std::numeric_limits<uint32_t>::max(), len=0;
            bool ok = get_col_span(std::string_view(gwas_lines.back()), idx_SNP, st, len);
            if (!ok) st = std::numeric_limits<uint32_t>::max();
            snp_span.emplace_back(st, len);
        }

        // 快解析 gCHR/gPOS/gA1/gA2 构建 gwas_vec（零拷贝 string_view
        std::string_view lv(gwas_lines.back());
        std::string_view vCHR, vPOS, vA1, vA2, vDummy;
        TabState st{0,0};

        int stop = std::max(std::max(gCHR, gPOS), std::max(gA1, gA2));
        scan_upto_col(lv, stop, gCHR, gPOS, gA1, gA2, -1, vCHR, vPOS, vA1, vA2, vDummy, st);

        int chr = canonical_chr_code_sv(vCHR);
        if (chr < 0) continue;

        int64_t pos = 0;
        if (!parse_i64(trim_ws(vPOS), pos) || pos <= 0) continue;

        AlleleKey ak = make_allele_key(trim_ws(vA1), trim_ws(vA2));
        if (ak.type == 2) continue;

        GWASRecord rec;
        rec.index  = idx;
        rec.chr    = chr;
        rec.pos    = pos;
        rec.allele = ak;

        gwas_vec.push_back(rec);
    }
    
    size_t n = gwas_lines.size();
    LOG_INFO("Loaded GWAS lines (data): " + std::to_string(n));

    //================ 基础 QC：过滤无效 N/beta/se/freq/P =================
    double maf   = P.maf_threshold;

    bool can_qc =  (idx_beta >= 0 ||
                    idx_se   >= 0 ||
                    idx_freq >= 0 ||
                    idx_pv   >= 0 ||
                    idx_n    >= 0);


    // 内部匹配循环用 uint8_t 更快；但 QC/dup 可能还用 vector<bool> -> 做一次性转换
    std::vector<bool> keep_qc_bool(n, true);
    std::vector<uint8_t> keep_qc_u8(n, 1);

    if (can_qc) {
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(gwas_lines, header, idx_beta, idx_se, idx_freq, idx_pv, idx_n, keep_qc_bool, maf);
        for (size_t i=0; i<n; ++i) keep_qc_u8[i] = keep_qc_bool[i] ? 1 : 0;
    } else {
        LOG_WARN("Cannot perform full QC in rsidImpu (missing beta/se/freq/N/p columns).");
    }

    //================ 构建排序用的 GWASRecord 向量 =================
    // 按 chr, pos 排序（稳定排序）
    std::stable_sort(gwas_vec.begin(), gwas_vec.end(),
        [](const GWASRecord& a, const GWASRecord& b){
            if (a.chr != b.chr) return a.chr < b.chr;
            return a.pos < b.pos;
        }
    );

    LOG_INFO("GWAS records sorted by CHR:POS for two-pointer matching.");

    //================ 准备匹配结果容器 =================
    // vector<uint8_t> 替代 vector<bool>
    std::vector<uint8_t> keep_u8(n, 0); // 是否最终进入主输出
    std::vector<std::string> rsid_vec(n); // 匹配到的 rsID

    //================ 单通扫描 dbSNP（Two-pointer merge） =================
    LineReader dbr(P.dbsnp_file);
    std::string dline;

    bool is_bim = ends_with(P.dbsnp_file, ".bim") || ends_with(P.dbsnp_file, ".bim.gz");

    int dCHR, dPOS, dA1, dA2, dRS;
    std::vector<std::string> dhdr;

    if (!is_bim) {
        // 有 header 的一般表格格式
        if (!dbr.getline(dline)) {
            LOG_ERROR("Empty dbSNP file.");
            exit(1);
        }
        strip_cr_inplace(dline);

        dhdr = split_tab(dline);
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

    // 真实扫描行数（而不是“命中候选”的行数）
    uint64_t scanned_total = 0;
    uint64_t scanned_valid_chrpos = 0;

    // 两段式解析：先只解析 CHR/POS，不命中就不解析 A1/A2/RS
    int stop_min = std::max(dCHR, dPOS);
    int stop_all = std::max({dCHR, dPOS, dA1, dA2, dRS});

    while (dbr.getline(dline)){
        if (dline.empty()) continue;
        strip_cr_inplace(dline);
        ++scanned_total;

        std::string_view lv(dline);

        std::string_view vCHR, vPOS, vA1, vA2, vRS;
        TabState st{0,0};

        // 1) 只扫到 CHR/POS（最大列号 stop_min）
        scan_upto_col(lv, stop_min, dCHR, dPOS, dA1, dA2, dRS, vCHR, vPOS, vA1, vA2, vRS, st);

        int dchr = canonical_chr_code_sv(vCHR);
        if (dchr < 0) continue;

        int64_t dpos = 0;
        if (!parse_i64(trim_ws(vPOS), dpos) || dpos <= 0) continue;

        ++scanned_valid_chrpos;

        // two-pointer 推进

        while (gi < Gn &&
            (gwas_vec[gi].chr < dchr ||
            (gwas_vec[gi].chr == dchr && gwas_vec[gi].pos < dpos))) {
            ++gi;
        }

        // 现在有两种可能：
        // 1) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos == dpos → 候选匹配
        // 2) gwas_vec[gi].chr == dchr && gwas_vec[gi].pos > dpos → 说明 dbSNP 这行在 GWAS 中不存在该 pos，继续读下一行 dbSNP
        if (gi >= Gn) break;

        if (gwas_vec[gi].chr != dchr || gwas_vec[gi].pos != dpos){
            // 2) 不命中：直接下一行（不解析 A1/A2/RS）
            if (scanned_total % 1000000ULL == 0) {
                LOG_INFO("[dbSNP two-pointer] scanned " + std::to_string(scanned_total/1000000ULL) + "M lines.");
            }
            continue;
        }

        // 3) 命中候选：再解析剩余列（A1/A2/RS）
        if (stop_all > stop_min){
            scan_upto_col(lv, stop_all, dCHR, dPOS, dA1, dA2, dRS, vCHR, vPOS, vA1, vA2, vRS, st);
        }

        // 等位基因规范化
        AlleleKey db_allele = make_allele_key(trim_ws(vA1), trim_ws(vA2));
        if (db_allele.type == 2) continue;

        // 可能有多个 GWAS 行在同一 chr:pos（或者多个 dbSNP 行同一 chr:pos）
        // 对所有该位置的 GWAS 进行尝试匹配
        size_t gj = gi;
        while (gj < Gn && 
            gwas_vec[gj].chr == dchr &&
            gwas_vec[gj].pos == dpos) {

            size_t orig_idx = gwas_vec[gj].index;
            // QC 未通过的行不做 rsID 匹配，但仍保留为 keep=false（之后输出到 unmatch）
            if (keep_qc_u8[orig_idx] &&
                gwas_vec[gj].allele.type == db_allele.type &&
                gwas_vec[gj].allele.key  == db_allele.key) {

                // 正向匹配 || 反向匹配
                keep_u8[orig_idx]     = 1;
                // rsID 只在命中时拷贝
                rsid_vec[orig_idx].assign(vRS.data(), vRS.size());
            }
            ++gj;
        }

        if(scanned_total % 1000000ULL == 0){
            LOG_INFO("[dbSNP two-pointer] scanned " + std::to_string(scanned_total/1000000ULL) + "M lines.");
        }
    }

    LOG_INFO("Two-pointer merge finished. dbSNP lines scanned: " + std::to_string(scanned_total) +
            ", valid CHR/POS lines: " + std::to_string(scanned_valid_chrpos));

    //================ 去重（按 rsID / P 值） =================
    if (P.remove_dup_snp) {
        std::vector<bool> keep_bool(n, false);
        for(size_t i=0; i<n; ++i) keep_bool[i] = (keep_u8[i] != 0);

        gwas_remove_dup(
            gwas_lines,
            header,
            idx_pv,
            rsid_vec,
            keep_bool
        );
        for (size_t i=0;i<n;++i) keep_u8[i] = keep_bool[i] ? 1 : 0;
    }
    
    //================ Writer（自动 txt / gz） =================
    bool out_is_gz = ends_with(P.out_file, ".gz");

    std::string out_main    = P.out_file;
    std::string base        = P.out_file;
    if (out_is_gz && base.size() > 3) {
        base.erase(base.size() - 3);  // 去掉 ".gz"
    }
    std::string out_unmatch = out_is_gz ? (base + ".unmatch.gz") : (P.out_file + ".unmatch");
    
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
        std::string h;
        for (size_t j=0; j<header.size(); j++) {
            if (j) h += "\t";
            h += header[j];
        }
        if (!has_SNP) h += "\tSNP";
        fout.write_line(h);
    } else {
        // 格式化 header
        std::string h;
        for (size_t j=0; j<spec.cols.size(); j++) {
            if (j) h += "\t";
            h += spec.cols[j];
        }
        fout.write_line(h);
    }

    // ================= 9) 输出 =================
    // 不再每行建 unordered_map；使用 FormatEngine fast path

    // 计算扫描 stop_col（只在 format != gwas 时用）
    int stop_out = gA1;
    stop_out = std::max(stop_out, gA2);
    stop_out = std::max(stop_out, idx_freq);
    stop_out = std::max(stop_out, idx_beta);
    stop_out = std::max(stop_out, idx_se);
    stop_out = std::max(stop_out, idx_pv);
    stop_out = std::max(stop_out, idx_n);
    
    for(size_t i=0; i<n; i++){
        if (!keep_u8[i]){
            funm.write_line(gwas_lines[i]);
            continue;
        }

        if (P.format == "gwas"){
            if (has_SNP) {
                // 直接用 span 替换 SNP 列（避免每次 find tab）
                if (P.format == "gwas" && has_SNP && i < snp_span.size()) {
                    auto [st, len] = snp_span[i];
                    if (st != std::numeric_limits<uint32_t>::max()) {
                        gwas_lines[i].replace((size_t)st, (size_t)len, rsid_vec[i]);
                    } else {
                        // ✅ 兜底：再算一次 span（防止预计算失败）
                        uint32_t st2=0, len2=0;
                        if (get_col_span(std::string_view(gwas_lines[i]), idx_SNP, st2, len2)) {
                            gwas_lines[i].replace((size_t)st2, (size_t)len2, rsid_vec[i]);
                        } else {
                            // ✅ 最终兜底：慢一点但不会错
                            replace_nth_column_inplace(gwas_lines[i], idx_SNP, rsid_vec[i]);
                        }
                    }
                }
                fout.write_line(gwas_lines[i]);
            } else {
                fout.write_line(gwas_lines[i] + "\t" + rsid_vec[i]);
            }
            continue;
        }

        // format != gwas: 快解析需要的列
        std::string_view lv(gwas_lines[i]);
        std::string_view vA1, vA2, vFreq, vBeta, vSe, vP, vN;
        TabState st{0,0};

        scan_upto_col7(lv, stop_out,
                       gA1, gA2, idx_freq, idx_beta, idx_se, idx_pv, idx_n,
                       vA1, vA2, vFreq, vBeta, vSe, vP, vN, st);

        // 使用 FormatEngine fast path
        FormatEngine::RowView row;
        row.SNP  = {std::string_view(rsid_vec[i]), true};
        row.A1   = {trim_ws(vA1), true};
        row.A2   = {trim_ws(vA2), true};

        row.freq = {trim_ws(vFreq), idx_freq >= 0};
        row.beta = {trim_ws(vBeta), idx_beta >= 0};
        row.se   = {trim_ws(vSe),   idx_se   >= 0};
        row.p    = {trim_ws(vP),    idx_pv   >= 0};
        row.N    = {trim_ws(vN),    idx_n    >= 0};
        
        fout.write_line(FE.format_line_fast(spec, row)); 
    }
}