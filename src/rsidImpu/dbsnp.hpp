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

using ChrMap = std::unordered_map<std::string, std::string>;
/**
 * Streaming-load dbSNP: only load the specified chromosome.
 * 
 * @param P     Args_RsidImpu parameters
 * @param CHR   canonical chromosome name ("1","2","X","Y","MT")
 * @return      unordered_map<key -> rsid>
 */
ChrMap streaming_load_chr(const Args_RsidImpu& P, const std::string& CHR);

#endif