#include "utils/FormatEngine.hpp"
#include <sstream>
#include <stdexcept>
FormatEngine::FormatEngine() {
    // ------- 原始 GWAS -------
    {
        FormatSpec spec;
        spec.name          = "gwas";
        spec.cols          = {};          // 由调用方自己决定如何输出（通常使用原始header + SNP）
        spec.required_rsid = false;
        spec.required_beta = false;
        spec.required_se   = false;
        spec.required_freq = false;
        spec.required_N    = false;
        spec.allow_missing = true;
        formats[spec.name] = spec;
    }

    // ------- COJO -------
    {
        FormatSpec spec;
        spec.name          = "cojo";
        spec.cols          = {"SNP","A1","A2","freq","b","se","p","N"};
        spec.required_rsid = true;   // 需要SNP
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_freq = true;
        spec.required_N    = true;
        spec.allow_missing = false;
        formats[spec.name] = spec;
    }

    // ------- Popcorn -------
    {
        FormatSpec spec;
        spec.name          = "popcorn";
        // 这里先放一个示例列集合，你后续可以按照Popcorn的文档改
        spec.cols          = {"SNP","A1","A2","beta","se","N"};
        spec.required_rsid = true;
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_N    = true;
        spec.required_freq = false;
        spec.allow_missing = false;
        formats[spec.name] = spec;
    }
    
    // ------- MR-MEGA -------
    {
        FormatSpec spec;
        spec.name          = "mrmega";
        // MR-MEGA一般要求beta,se,z,N,p等等，这里先给一个基础示例
        spec.cols          = {"SNP","A1","A2","beta","se","P","N"};
        spec.required_rsid = true;
        spec.required_beta = true;
        spec.required_se   = true;
        spec.required_N    = true;
        spec.required_freq = false;
        spec.allow_missing = false;
        formats[spec.name] = spec;
    }
}

FormatSpec FormatEngine::get_format(const std::string& name) const{
    auto it = formats.find(name);
    if (it == formats.end()){
        throw std::runtime_error("Unknown GWAS format: " + name);
    }

    return it->second;
}

std::string FormatEngine::format_line(
    const FormatSpec& spec,
    const std::unordered_map<std::string,std::string>&row) const
{
    std::ostringstream oss;
    for (size_t i = 0; i < spec.cols.size(); i++){
        const auto &col = spec.cols[i];
        if (i) oss << "\t";

        auto it = row.find(col);
        if (it == row.end()) {
            if (spec.allow_missing) {
                oss << "";
            } else {
                // throw std::runtime_error("Missing required column in row: " + col);
                oss << "";
            }
        } else {
            oss << it->second;
        }
    }
    return oss.str();
}