//
//  allele.cpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/util.hpp"
#include "rsidImpu/allele.hpp"

using namespace std;

// A↔T / C↔G 互补
char comp_base(char b){
    switch (toupper((unsigned char)b)){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
    }
    return b;
}

// 多字符等位基因（如 indel）互补
string comp_allele(const string &a){
    string t = a;
    for (auto &c : t) c = comp_base(c);
    return t;
}

pair<string,string> canonical_alleles(const string &a1, const string &a2){
    string n1 = upper(trim(a1));
    string n2 = upper(trim(a2));

    string f1=n1, f2=n2;
    if (f2 < f1) swap(f1, f2);

    string c1=comp_allele(n1), c2=comp_allele(n2);
    string cf1=c1, cf2=c2;
    if (cf2 < cf1) swap(cf1, cf2);

    return (cf1 + ":" + cf2 < f1 + ":" + f2) ? make_pair(cf1,cf2)
                                            : make_pair(f1,f2);
}

string make_key(const string &chr, const string &pos, const string &a1, const string &a2){
    auto canon = canonical_alleles(a1, a2);
    return trim(chr) + ":" + trim(pos) + ":" + canon.first + ":" + canon.second;
}