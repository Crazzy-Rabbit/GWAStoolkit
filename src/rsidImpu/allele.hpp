//
//  rsidImpu.cpp
//
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_ALLELE_HPP
#define RSIDIMPU_ALLELE_HPP

#include <string>
#include <utility>

char comp_base(char b);
std::string comp_allele(const std::string&);
std::pair<std::string,std::string> canonical_alleles(const std::string&, const std::string&);
std::string make_key(const std::string&, const std::string&, const std::string&, const std::string&);

#endif
