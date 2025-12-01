#ifndef GWASTOOLKIT_FORMAT_ENGINE_HPP
#define GWASTOOLKIT_FORMAT_ENGINE_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <map>

struct FormatSpec {
    std::string name;
    std::vector<std::string> cols;
    bool required_rsid = false;
    bool required_beta = false;
    bool required_se   = false;
    bool required_freq = false;
    bool required_N    = false;
    bool allow_missing = false;
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

private:
    std::map<std::string,FormatSpec> formats;

    // 核心：将各种格式的列名映射到内部 key
    // "b" → "beta"     "BETA" → "beta"
    // "FREQ" → "freq"  "SE" → "se"
    std::string normalize_key(const std::string& col) const;
    static std::string to_lower(std::string s);
};

#endif