//
//  args.hpp
//  rsidImpu
//
//  Created by Lulu Shi on 28/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//

#ifndef RSIDIMPU_ARGS_HPP
#define RSIDIMPU_ARGS_HPP

#include <string>
#include <vector>
#include <map>

// ----------------------【公共字段】-------------------------
struct CommonArgs {
    std::string gwas_file;
    std::string out_file;

        // ---- 通用 GWAS 列名参数 ----
    std::string col_SNP  = "SNP";
    std::string g_chr    = "CHR";   // --chr
    std::string g_pos    = "POS";   // --pos
    std::string g_A1     = "A1";    // --A1
    std::string g_A2     = "A2";    // --A2
    std::string g_p      = "p";     // --pval
    std::string col_freq = "freq";
    std::string col_beta = "b";
    std::string col_se   = "se";
    std::string col_n    = "N";

    std::string format   = "gwas";

    bool remove_dup_snp  = false;
    double maf_threshold = 0.01;

    int threads          = 1;
    bool log_enabled     = false;
    std::string log_file;
};

// ----------------------【rsid-impu 子命令专用】-------------------------
struct Args_RsidImpu : public CommonArgs {
    std::string dbsnp_file;
    std::string d_chr;
    std::string d_pos;
    std::string d_A1;
    std::string d_A2;
    std::string d_rsid;
};

// ----------------------【convert 子命令专用】-------------------------
struct Args_Convert : public CommonArgs {
    //
};

// ----------------------【or2beta 子命令专用】-------------------------
struct Args_Or2Beta : public CommonArgs {
    std::string col_or = "OR";
};

struct Args_CalNeff : public CommonArgs {
    bool is_single = false;
    bool is_column = false;

    int case_n = 0;
    int control_n = 0;

    std::string case_col;
    std::string control_col;
};

// ----------------------【解析器接口】-------------------------
void print_rsidimpu_help();
void print_convert_help();
void print_or2beta_help();
void print_calneff_help();

Args_RsidImpu  parse_args_rsidimpu(int argc, char* argv[]);
Args_Convert   parse_args_convert(int argc, char* argv[]);
Args_Or2Beta   parse_args_or2beta(int argc, char* argv[]);
Args_CalNeff  parse_args_calneff(int argc, char* argv[]);

#endif
