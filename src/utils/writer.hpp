//
//  writer.hpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef TOOLKIT_WRITER_HPP
#define TOOLKIT_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <zlib.h>

// 可指定输出格式
// 如果文件名以 ".gz" 结尾 → 用 gzopen 写 gzip
// 否则 → 用 ofstream 写普通文本

class Writer {
public:
    Writer(const std::string &filename, const std::string &format = "gwas");
    ~Writer();

    void write_line(const std::string &line);
    bool good() const { return ok_; }

private:
    bool use_gz_ = false;
    bool ok_ = false;

    std::ofstream ofs_;
    gzFile gzfp_ = nullptr;
};

#endif