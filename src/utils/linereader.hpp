//
//  linereader.hpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_LINEREADER_HPP
#define RSIDIMPU_LINEREADER_HPP

#include <string>

class LineReader {
public:
    LineReader(const std::string&);
    ~LineReader();
    bool getline(std::string &line);

private:
    std::string fname;
    bool gz;
    void* gzfp;
    std::ifstream* fin;

    static bool ends_with(const std::string&, const std::string&);
};

#endif
