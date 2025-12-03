//
//  rsidImpu.cpp
//
//  Created by Lulu Shi on 24/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#include "utils/util.hpp"
#include "utils/log.hpp"

#include <sstream>
#include <algorithm>
#include <cctype>

using namespace std;

vector<string> split(const string &s){
    vector<string> out;
    string token;
    stringstream ss(s);
    while (ss >> token) out.push_back(token);
    return out;
}

/***********************
 * trim(): remove spaces + tabs + \r + \n
 ***********************/
string trim(const string &s){
    size_t i = 0, j = s.size();
    while (i < j && (isspace((unsigned char)s[i]) || s[i] == '\r' || s[i] == '\n'))
        i++;
    while (j > i && (isspace((unsigned char)s[j - 1]) || s[j - 1] == '\r' || s[j - 1] == '\n'))
        j--;
    return s.substr(i, j - i);
}

// 转成大写
string upper(const string &s){
    string t = s;
    for (auto &c : t) c = toupper((unsigned char)c);
    return t;
}

string norm_chr(const string &chr){
    string c = trim(chr);

    // 去掉 "chr" 前缀 (Chr / CHR / chr)
    if (c.size() >= 3) {
        string prefix = c.substr(0, 3);
        transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
        if (prefix == "chr"){
            return c.substr(3);
        }
    }
    return c;
}

// 后缀检查
bool ends_with(const string &s, const string &suffix){
    if (s.size() < suffix.size()) return false;
    return equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

// 
int find_col(const vector<string> &header, const string &colname){
    for (int i = 0; i < (int)header.size(); i++)
        if (header[i] == colname) return i;
    return -1;
}

//
void require(bool cond, const std::string& msg){
    if(!cond){
        LOG_ERROR(msg);
        exit(1);
    }
}


static inline bool starts_with(const std::string& s, const std::string& p){
    return s.rfind(p, 0) == 0;
}

std::string canonical_chr(const std::string& raw){
    std::string x = raw;
    
    // 去除回车键
    x.erase(std::remove(x.begin(), x.end(), '\r'), x.end());
    
    // 同意大小写
    for (auto& c : x) c = std::toupper(c);

    // 去掉 CHR前缀
    if (starts_with(x, "CHR"))
        x = x.substr(3);

    // RefSeq: NC_000001.11 -> 1
    if (starts_with(x, "NC_")){
        if (x.size() >= 9){
            std::string num = x.substr(3, 6);

            while (num.size() > 1 && num[0] == '0')
                num.erase(0,1);

            int chr = atoi(num.c_str());
            if(1 <= chr && chr <= 22)
                return std::to_string(chr);
            if (chr == 23)
                return "X";
            if (chr == 24)
                return "Y";
            if (chr == 12920)
                return "MT";
        }
    }
    
    // 去前导 0：01 → 1 
    while (x.size() > 1 && x[0] == '0')
        x.erase(0,1);

    // 常用别名
    if (x == "M" || x == "MT" || x == "MTDNA")
        return "MT";

    if (x == "23") return "X";
    if (x == "24") return "Y";

    return x;
}