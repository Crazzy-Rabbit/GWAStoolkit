#include "utils/FormatEngine.hpp"
#include <stdexcept>
#include <algorithm>

// 小工具：ASCII lower
static inline char low(char c){ return (char)std::tolower((unsigned char)c); }

FieldId FormatEngine::col_to_field_id(std::string_view col)
{
    // case-insensitive match for known fields
    // 允许 "b"/"beta", "se"/"SE", "freq"/"FREQ", "n"/"N"
    auto eq_ci = [](std::string_view a, std::string_view b){
        if (a.size() != b.size()) return false;
        for (size_t i=0;i<a.size();++i) if (low(a[i]) != low(b[i])) return false;
        return true;
    };

    if (eq_ci(col, "SNP"))  return FieldId::SNP;
    if (eq_ci(col, "A1"))   return FieldId::A1;
    if (eq_ci(col, "A2"))   return FieldId::A2;

    if (eq_ci(col, "freq")) return FieldId::FREQ;
    if (eq_ci(col, "b") || eq_ci(col, "beta") || eq_ci(col, "BETA")) return FieldId::BETA;
    if (eq_ci(col, "se") || eq_ci(col, "SE")) return FieldId::SE;
    if (eq_ci(col, "p") || eq_ci(col, "P")) return FieldId::P;
    if (eq_ci(col, "n") || eq_ci(col, "N")) return FieldId::N;

    return FieldId::UNKNOWN;
}

std::string_view FormatEngine::key_of(FieldId id)
{
    // row 内部 key（你原先的 normalize_key 目标）
    switch(id){
        case FieldId::SNP:  return "SNP";
        case FieldId::A1:   return "A1";
        case FieldId::A2:   return "A2";
        case FieldId::FREQ: return "freq";
        case FieldId::BETA: return "beta";
        case FieldId::SE:   return "se";
        case FieldId::P:    return "p";
        case FieldId::N:    return "N";
        default:            return "";
    }
}

// string append 替代 ostringstream
inline void FormatEngine::append_tabbed(std::string& out, bool& first, std::string_view sv)
{
    if (!first) out.push_back('\t');
    first = false;
    out.append(sv.data(), sv.size());
}

FormatEngine::FormatEngine() {
    // ------- 原始 GWAS -------
    {
        FormatSpec spec;
        spec.name          = "gwas";
        spec.cols          = {};          // 由调用方自己决定如何输出（通常使用原始header + SNP）
        spec.allow_missing = true;
        spec.field_ids     = {};
        formats[spec.name] = spec;
    }

    // ------- COJO -------
    {
        FormatSpec spec;
        spec.name          = "cojo";
        spec.cols          = {"SNP","A1","A2","freq","b","se","p","N"};
        spec.required_rsid = true;
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_freq = true;
        spec.required_N    = true;
        spec.allow_missing = false;

        // [OPT-FE-1] 预计算 field_ids
        spec.field_ids.reserve(spec.cols.size());
        for (auto &c : spec.cols) spec.field_ids.push_back(col_to_field_id(c));

        formats[spec.name] = spec;
    }

    // ------- Popcorn -------
    {
        FormatSpec spec;
        spec.name          = "popcorn";
        spec.cols          = {"SNP","A1","A2","freq","beta","SE","N"};
        spec.required_rsid = true;
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_freq = true;
        spec.required_N    = true;
        spec.allow_missing = false;

        spec.field_ids.reserve(spec.cols.size());
        for (auto &c : spec.cols) spec.field_ids.push_back(col_to_field_id(c));

        formats[spec.name] = spec;
    }
    
    // ------- MR-MEGA -------
    {
        FormatSpec spec;
        spec.name          = "mrmega";
        spec.cols          = {"SNP","A1","A2","FREQ","BETA","SE","P","N"};
        spec.required_rsid = true;
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_freq = true;
        spec.required_N    = true;
        spec.allow_missing = false;

        spec.field_ids.reserve(spec.cols.size());
        for (auto &c : spec.cols) spec.field_ids.push_back(col_to_field_id(c));

        formats[spec.name] = spec;
    }
}

FormatSpec FormatEngine::get_format(const std::string& name) const{
    auto it = formats.find(name);
    if (it == formats.end()){
        throw std::runtime_error("Unknown GWAS format: " + name);
    }
    return it->second; // 复制一份 spec（包含 field_ids）
}

std::string FormatEngine::format_line(
    const FormatSpec& spec,
    const std::unordered_map<std::string,std::string>&row) const
{
    // 不用 ostringstream
    std::string out;
    out.reserve(128);

    bool first = true;
    for (size_t i=0; i<spec.cols.size(); ++i){
        FieldId id = (i < spec.field_ids.size() ? spec.field_ids[i] : FieldId::UNKNOWN);

        std::string_view key = key_of(id);
        if (key.empty()) {
            // fallback: 找原列名（极少见）
            key = spec.cols[i];
        }

        auto it = row.find(std::string(key));
        if (it != row.end()) {
            append_tabbed(out, first, it->second);
        } else {
            if (!spec.allow_missing) {
                throw std::runtime_error("Missing required column [" + spec.cols[i] + "]");
            }
            append_tabbed(out, first, "");
        }
    }
    return out;
}

// ---------------- 新增 fast path：不走 unordered_map ----------------
std::string FormatEngine::format_line_fast(const FormatSpec& spec, const RowView& row) const
{
    std::string out;
    out.reserve(128);
    bool first = true;

    for (size_t i=0; i<spec.cols.size(); ++i){
        FieldId id = spec.field_ids[i];

        const CellView* cell = nullptr;
        switch(id){
            case FieldId::SNP:  cell = &row.SNP; break;
            case FieldId::A1:   cell = &row.A1;  break;
            case FieldId::A2:   cell = &row.A2;  break;
            case FieldId::FREQ: cell = &row.freq;break;
            case FieldId::BETA: cell = &row.beta;break;
            case FieldId::SE:   cell = &row.se;  break;
            case FieldId::P:    cell = &row.p;   break;
            case FieldId::N:    cell = &row.N;   break;
            default:            cell = nullptr;  break;
        }

        if (cell && cell->present) {
            append_tabbed(out, first, cell->v);
        } else {
            if (!spec.allow_missing) {
                throw std::runtime_error("Missing required column [" + spec.cols[i] + "]");
            }
            append_tabbed(out, first, "");
        }
    }
    return out;
}