//
//  writer.cpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/writer.hpp"
#include "utils/util.hpp"   // 用里面的 ends_with
#include "utils/log.hpp"

#include <iostream>

Writer::Writer(const std::string &filename, const std::string & /*format*/)
{
    // 判断是否 .gz 结尾
    if (ends_with(filename, ".gz")) {
        use_gz_ = true;
        gzfp_ = gzopen(filename.c_str(), "wb");
        if (!gzfp_) {
            LOG_ERROR("Error: cannot open gzip file for writing: " + filename);
            ok_ = false;
            return;
        }
    } else {
        ofs_.open(filename);
        if (!ofs_) {
            LOG_ERROR("Error: cannot open text file for writing: " + filename);
            ok_ = false;
            return;
        }
    }

    ok_ = true;
}

Writer::~Writer()
{
    if (use_gz_) {
        if (gzfp_) gzclose(gzfp_);
    } else {
        if (ofs_.is_open()) ofs_.close();
    }
}

void Writer::write_line(const std::string &line)
{
    if (!ok_) return;

    if (use_gz_) {
        gzputs(gzfp_, (line + "\n").c_str());
    } else {
        ofs_ << line << '\n';
    }
}