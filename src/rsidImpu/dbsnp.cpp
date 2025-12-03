//
//  dbsnp.cpp
//  GWAStoolkit
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

/**
 * Streaming-load dbSNP by chromosome.
 * Only lines with CHR==target chromosome are loaded.
 */

ChrMap streaming_load_chr(const Args_RsidImpu& P, const std::string& CHR)
{
    ChrMap chrmap;
    chrmap.reserve(3000000);   // pre-allocate to speed up

    LineReader reader(P.dbsnp_file);
    string line;

    bool is_bim = ends_with(P.dbsnp_file,".bim") ||
                    ends_with(P.dbsnp_file,".bim.gz");
    
    int dCHR, dPOS, dA1, dA2, dRS;

    // parse header
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

    // streaming reading
    size_t count = 0;
    size_t matched = 0;

    LOG_INFO("[dbSNP] Loading chromosome " + CHR + " ...");

    while (reader.getline(line)) {
        if (line.empty()) continue;
        count++;

        if (count % 1000000 == 0){
            LOG_INFO("[dbSNP:" + CHR + "] scanned " + to_string(count/1000000) + "M lines");
        }

        auto f = split(line);

        // string chr = norm_chr(f[dCHR]);
        string chr = canonical_chr(f[dCHR]);
        if (chr != CHR) continue;   // <--- 仅加载该chr的记录

        matched++;

        string pos = f[dPOS];
        string a1  = f[dA1];
        string a2  = f[dA2];
        string rs  = f[dRS];

        auto canon = canonical_alleles(a1,a2);
        string key = pos + ":" + canon.first + ":" + canon.second;

        chrmap[key] = rs;
    }
    LOG_INFO("[dbSNP:" + CHR + "] finish. scanned " + to_string(count) 
            + ", matched " + to_string(matched));
    return chrmap;
}