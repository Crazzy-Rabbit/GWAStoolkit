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
#include <cctype>   // [LOWERCASE] for std::toupper

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

// [LOWERCASE] 允许小写输入（很多 GWAS/db 文件会有小写）
inline char to_upper_base(char c){
    return (char)std::toupper((unsigned char)c);
}


inline bool is_acgt(char c){
    c = to_upper_base(c);   
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

// pack sorted codes
inline uint64_t pack_snp_pair(uint8_t x, uint8_t y){
    if (x > y) std::swap(x, y);
    return (uint64_t(x) << 4) | y;
}

// [STRAND-FLIP] 互补链不变的 SNP key：min(原始, 互补)
inline uint64_t snp_pair_strand_invariant(char a1, char a2){
    uint8_t x = snp_code(a1);
    uint8_t y = snp_code(a2);
    // 调用前通常已 is_acgt 检查，这里再防御一下
    if (x > 3 || y > 3) return 0;

    uint64_t k1 = pack_snp_pair(x, y);

    // A(0)<->T(3), C(1)<->G(2) 正好是 code' = 3 - code
    uint64_t k2 = pack_snp_pair(uint8_t(3 - x), uint8_t(3 - y));

    return std::min(k1, k2);
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

        // [STRAND-FLIP] 这里改用 strand-invariant key
        return {0, snp_pair_strand_invariant(a1[0], a2[0])};
    }

    // INDEL / multi-base ACGT
    // 注：这里暂不做反向互补（indel 表达复杂且会影响性能/一致性）
    if (is_acgt_string(a1) && is_acgt_string(a2)) {
        return {1, indel_pair(a1, a2)};
    }

    // OTHER (rare, ignored in matching)
    return {2, 0};
}