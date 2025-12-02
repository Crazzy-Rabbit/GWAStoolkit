#include "utils/FormatEngine.hpp"
#include <sstream>
#include <stdexcept>
#include <algorithm>

// to lower
std::string FormatEngine::to_lower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                    [](unsigned char c){return std::tolower(c); });
    return s;
}
//
std::string FormatEngine::normalize_key(const std::string& col) const{
    std::string low = to_lower(col);

    if(low == "b" || low == "beta") return "beta";
    if(low == "se")                 return "se";
    if(low == "p")                  return "p";
    if(low == "freq")               return "freq";
    if(low == "n")                  return "N";

    if (col == "SNP" || col == "A1" || col == "A2")
        return col;

    return col;
}


FormatEngine::FormatEngine() {
    // ------- 原始 GWAS -------
    {
        FormatSpec spec;
        spec.name          = "gwas";
        spec.cols          = {};          // 由调用方自己决定如何输出（通常使用原始header + SNP）
        spec.allow_missing = true;
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
        if (i) oss << "\t";

        std::string col = spec.cols[i];

        // tranform to inhouse key
        std::string key = normalize_key(col);

        auto it = row.find(key);
        if (it != row.end()) {
            oss << it->second;
        } else  {
            if (!spec.allow_missing){
                throw std::runtime_error(
                    "Missing required column [" + col + "], internal key [" + key + "] not provided!"
                );
            }
            oss << "";
        }
    }

    return oss.str();
}