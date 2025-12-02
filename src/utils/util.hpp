//
//  util.hpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_UTIL_HPP
#define RSIDIMPU_UTIL_HPP

#include <string>
#include <vector>

std::vector<std::string> split(const std::string& s);
std::string trim(const std::string& s);
std::string upper(const std::string& s);
std::string norm_chr(const std::string &chr);
bool ends_with(const std::string& s, const std::string& suffix);
int find_col(const std::vector<std::string>& header, const std::string& colname);
void require(bool cond, const std::string& msg);
#endif
