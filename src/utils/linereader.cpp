//
//  linereader.cpp
//  rsidImpu
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "linereader.hpp"
#include <zlib.h>
#include <fstream>

using namespace std;

LineReader::LineReader(const string &filename){
    fname = filename;
    gz = false;
    gzfp = nullptr;
    fin = nullptr;

    if (ends_with(fname, ".gz")){
        gz = true;
        gzfp = gzopen(fname.c_str(), "rb");
        if (!gzfp) throw runtime_error("Cannot open gz file: " + fname);
    } else {
        fin = new ifstream(fname);
        if (!fin->good()){
            delete fin;
            fin = nullptr;
            throw runtime_error("Cannot open file: " + fname);
        }
    }
}

LineReader::~LineReader(){
    if (gz) {
        if (gzfp) gzclose((gzFile)gzfp);
    } else {
        if (fin) {
            fin->close();
            delete fin;
        }
    }
}

bool LineReader::getline(string &line){
    if (gz) {
        const int BUF_SIZE = 1<<16;
        static thread_local char buf[BUF_SIZE];
        char *res = gzgets((gzFile)gzfp, buf, BUF_SIZE);
        if (!res) return false;
        line.assign(buf);
        while (!line.empty() && (line.back()=='\n' || line.back()=='\r'))
            line.pop_back();
        return true;
    }
    return (bool)std::getline(*fin, line);
}

bool LineReader::ends_with(const string &s, const string &suffix){
    if (s.size() < suffix.size()) return false;
    return equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}