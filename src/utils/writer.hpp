//
//  writer.hpp
//  rsidImpu
//  Created by Lulu Shi on 25/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_WRITER_HPP
#define RSIDIMPU_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <zlib.h>

// 可指定输出格式
// 如果文件名以 ".gz" 结尾 → 用 gzopen 写 gzip
// 否则 → 用 ofstream 写普通文本

class Writer {
public:
    Writer(const std::string &filename, const std::string &format);
    ~Writer();

    void write_line(const std::string &line);

    void write_cojo_header();

    bool good() const { return ok_; }

private:
    bool use_gz_;
    bool ok_;
    std::string format_;  // "gwas", "smr", etc.
    std::ofstream ofs_;
    gzFile gzfp_;
};

#endif