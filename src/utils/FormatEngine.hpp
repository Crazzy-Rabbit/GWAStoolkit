#ifndef GWASTOOLKIT_FORMAT_ENGINE_HPP
#define GWASTOOLKIT_FORMAT_ENGINE_HPP

#include <string>
#include <vector>
#include <unordered_map>

enum class GWASFormat {
    GWAS,
    COJO,
    SMR,
    LDSC,
    Popcorn,
    MR_MEGA
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
};

class FormatEngine {
public:
    FormatEngine();
    FormatSpec get_format(const std::string& name) const;
    std::string format_line(const FormatSpec& spec,
                            const std::unordered_map<std::string,std::string>& row) const;

private:
    std::unordered_map<std::string,FormatSpec> formats;
};

#endif