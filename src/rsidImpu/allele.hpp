//
//  rsidImpu.cpp
//
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#pragma once
#include <string_view>
#include <cstdint>
#include <algorithm>

// =======================================================
// AlleleKey: compact allele representation
// type: 0=SNP, 1=INDEL, 2=OTHER
// =======================================================

struct AlleleKey {
    uint8_t  type;
    uint64_t key;

    inline bool operator==(const AlleleKey& o) const {
        return type == o.type && key == o.key;
    }
};

// ---------------- basic utils ----------------

inline bool is_acgt(char c){
    return c=='A' || c=='C' || c=='G' || c=='T';
}

inline bool is_acgt_string(std::string_view s){
    for (char c : s)
        if (!is_acgt(c)) return false;
    return true;
}

// ---------------- SNP fast path ----------------

inline uint8_t snp_code(char a){
    switch(a){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 255;
    }
}

inline uint64_t snp_pair(char a1, char a2){
    uint8_t x = snp_code(a1);
    uint8_t y = snp_code(a2);
    if (x > y) std::swap(x,y);
    return (uint64_t(x) << 4) | y;
}

// ---------------- INDEL / multi-base ----------------

inline uint64_t hash_sv(std::string_view s){
    uint64_t h = 1469598103934665603ULL; // FNV-1a
    for (char c : s){
        h ^= (uint64_t)c;
        h *= 1099511628211ULL;
    }
    return h;
}

inline uint64_t indel_pair(std::string_view a1, std::string_view a2){
    if (a1 > a2) std::swap(a1, a2);
    uint64_t h1 = hash_sv(a1);
    uint64_t h2 = hash_sv(a2);
    return h1 ^ (h2 << 1);
}

// ---------------- unified entry ----------------

inline AlleleKey make_allele_key(std::string_view a1,
                                std::string_view a2)
{
    // SNP
    if (a1.size()==1 && a2.size()==1 &&
        is_acgt(a1[0]) && is_acgt(a2[0])) {
        return {0, snp_pair(a1[0], a2[0])};
    }

    // INDEL / multi-base ACGT
    if (is_acgt_string(a1) && is_acgt_string(a2)) {
        return {1, indel_pair(a1, a2)};
    }

    // OTHER (rare, ignored in matching)
    return {2, 0};
}