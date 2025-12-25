#pragma once
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <cstdint>

enum class FieldId : uint8_t {
    SNP, A1, A2, FREQ, BETA, SE, P, N, UNKNOWN
};

struct FormatSpec {
    std::string name;
    std::vector<std::string> cols;
    bool required_rsid = false;
    bool required_beta = false;
    bool required_se   = false;
    bool required_freq = false;
    bool required_N    = false;
    bool allow_missing = false;

    std::vector<FieldId> field_ids;
};

class FormatEngine {
public:
    FormatEngine();
    // 获取指定格式
    FormatSpec get_format(const std::string& name) const;
    // 映射输出
    std::string format_line(
        const FormatSpec& spec,
        const std::unordered_map<std::string,std::string>& row
    ) const;
    // fast path
    struct CellView {
        std::string_view v;
        bool present = false;
    };
    struct RowView {
        CellView SNP, A1, A2, freq, beta, se, p, N;
    };

    std::string format_line_fast(const FormatSpec& spec, const RowView& row) const;

private:
    std::unordered_map<std::string, FormatSpec> formats;

    static FieldId col_to_field_id(std::string_view col);
    static std::string_view key_of(FieldId id);

    static inline void append_tabbed(std::string& out, bool& first, std::string_view sv);
};