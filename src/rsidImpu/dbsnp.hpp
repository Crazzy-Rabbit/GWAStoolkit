//
//  dbsnp.hpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//


#ifndef RSIDIMPU_DBSNP_HPP
#define RSIDIMPU_DBSNP_HPP

#include "utils/args.hpp"

#include <unordered_map>
#include <string>

// std::unordered_map<std::string,std::string>
using DBMap = std::unordered_map<
    std::string,
    std::unordered_map<std::string, std::string>
>;
DBMap load_dbsnp(const Args_RsidImpu& P);

#endif