//
//  dbsnp.cpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/linereader.hpp"
#include "rsidImpu/dbsnp.hpp"
#include "rsidImpu/allele.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

DBMap load_dbsnp(const Args_RsidImpu& P) {
    DBMap mapdb;   // chr → (key → rsid)

    LineReader reader(P.dbsnp_file);
    string line;

    bool is_bim = ends_with(P.dbsnp_file,".bim") ||
                    ends_with(P.dbsnp_file,".bim.gz");

    int dCHR, dPOS, dA1, dA2, dRS;

    if (!is_bim){
        reader.getline(line);
        auto hdr = split(line);
        dCHR = find_col(hdr, P.d_chr);
        dPOS = find_col(hdr, P.d_pos);
        dA1  = find_col(hdr, P.d_A1);
        dA2  = find_col(hdr, P.d_A2);
        dRS  = find_col(hdr, P.d_rsid);

        if (dCHR<0||dPOS<0||dA1<0||dA2<0||dRS<0){
            LOG_ERROR("dbSNP header incomplete.");
            exit(1);
        }
    } else {
        dCHR=0; dRS=1; dPOS=3; dA1=4; dA2=5;
    }

    while (reader.getline(line)) {
        if (line.empty()) continue;

        auto f = split(line);

        string chr = norm_chr(f[dCHR]);
        string pos = f[dPOS];
        string a1  = f[dA1];
        string a2  = f[dA2];
        string rs  = f[dRS];

        // 用精简的 key: pos + alleles
        auto canon = canonical_alleles(a1,a2);
        string key = pos + ":" + canon.first + ":" + canon.second;

        mapdb[chr][key] = rs;
    }

    LOG_INFO("Loaded dbSNP entries: " + std::to_string(mapdb.size()) + " chromosomes.");
    return mapdb;
}