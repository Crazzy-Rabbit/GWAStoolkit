//
//  args.cpp
//  rsidImpu
//
//  Created by Lulu Shi on 28/11/2025.
//  Copyright © 2025 Lulu Shi. All rights reserved.
//


#include "utils/args.hpp"
#include "utils/log.hpp"

#include <iostream>
#include <cstdlib>
#include <set>

using namespace std;

static const set<string> supported_formats = {
    "gwas", "cojo", "popcorn", "mrmega"
};

// ★ 所有合法的参数名（包括 flag 和带值项）
static const set<string> common_params = {
    "--gwas-summary", "--out",
    "--SNP", "--chr", "--pos", "--A1", "--A2", "--pval",
    "--freq", "--beta", "--se", "--n",
    "--format",
    "--maf", "--remove-dup-snp",
    "--threads", "--log"
};
static const std::set<std::string> rsidimpu_params = {
    "--dbsnp", "--dbchr", "--dbpos", "--dbA1", "--dbA2", "--dbrsid"
};
static const std::set<std::string> convert_params = {};
static const std::set<std::string> or2beta_params = {
    "--or"
};

// =============== 通用错误检查 ===================
static void require(bool cond, const string& msg){
    if (!cond) {
        LOG_ERROR(msg);
        exit(1);
    }
}

static void check_format(const string& fmt){
    if (!supported_formats.count(fmt)) {
        LOG_ERROR("Unsupported format: " + fmt + 
        " (supported: gwas, cojo, popcorn, mrmega)");
        exit(1);
    }
}

// ------------------------- 公共解析  ---------------------
static void parse_common(CommonArgs& C, map<string,string>& args) {
    require(args.count("--gwas-summary"), "Missing required: --gwas-summary");
    require(args.count("--out"), "Missing required: --out");

    C.gwas_file = args["--gwas-summary"];
    C.out_file  = args["--out"];

    if (args.count("--threads"))
        C.threads = stoi(args["--threads"]);

    if (args.count("--log")) {
        C.log_enabled = true;
        C.log_file = args["--log"];
    }

    if (args.count("--remove-dup-snp")){
        C.remove_dup_snp = true;
    }
    
    if (args.count("--maf")) {
        C.maf_threshold = stod(args["--maf"]);
    }

    // 公共列名
    if (args.count("--SNP"))   C.col_SNP  = args["--SNP"];
    if (args.count("--chr"))   C.g_chr    = args["--chr"];
    if (args.count("--pos"))   C.g_pos    = args["--pos"];
    if (args.count("--A1"))    C.g_A1     = args["--A1"];
    if (args.count("--A2"))    C.g_A2     = args["--A2"];
    if (args.count("--pval"))  C.g_p      = args["--pval"];
    if (args.count("--freq"))  C.col_freq = args["--freq"];
    if (args.count("--beta"))  C.col_beta = args["--beta"];
    if (args.count("--se"))    C.col_se   = args["--se"];
    if (args.count("--n"))     C.col_n    = args["--n"];

    if (args.count("--format"))
        C.format = args["--format"];
    check_format(C.format);
}

// ======================================================
//                     HELP 信息
// ======================================================
void print_rsidimpu_help() {
    cerr <<
    "Usage: GWAStoolkit rsidImpu [options]\n\n"
    "Required:\n"
    "  --gwas-summary FILE        Input GWAS summary statistics\n"
    "  --dbsnp FILE               dbSNP or BIM file\n"
    "  --out FILE                 Output file\n"
    "  --dbchr COL --dbpos COL --dbA1 COL --dbA2 COL --dbrsid COL\n\n"
    "Required GWAS columns:\n"
    "  --chr COL   (default: CHR)\n"
    "  --pos COL   (default: POS)\n"
    "  --A1  COL   (default: A1)\n"
    "  --A2  COL   (default: A2)\n\n"
    "Optional format:\n"
    "  --format gwas|cojo|popcorn|mrmega   Output format (default: gwas)\n\n"
    "Optional columns:\n"
    "  --pval COL   (default: p)\n"
    "  --freq COL   (default: freq)\n"
    "  --beta COL   (default: b)\n"
    "  --se   COL   (default: se)\n"
    "  --n    COL   (default: N)\n\n"
    "QC options:\n"
    "  --maf VAL              MAF threshold (default: 0.01)\n"
    "  --remove-dup-snp       Keep smallest-P SNP when duplicates exist\n\n"
    "Other options:\n"
    "  --threads N\n"
    "  --log FILE\n";
}

void print_convert_help() {
    cerr <<
    "Usage: GWAStoolkit convert [options]\n\n"
    "Required:\n"
    "  --gwas-summary FILE\n"
    "  --out FILE\n"
    "  --format gwas|cojo|popcorn|mrmega\n"
    "  --SNP COL        SNP identifier column\n\n"
    "Required GWAS columns for conversion:\n"
    "  --A1 COL   (default: A1)\n"
    "  --A2 COL   (default: A2)\n"
    "  --pval COL (default: p)\n"
    "  --freq COL (default: freq)\n"
    "  --beta COL (default: b)\n"
    "  --se   COL (default: se)\n"
    "  --n    COL (default: N)\n\n"
    "QC options:\n"
    "  --maf VAL\n"
    "  --remove-dup-snp\n\n"
    "Other options:\n"
    "  --threads N\n"
    "  --log FILE\n";
}

void print_or2beta_help() {
    cerr <<
    "Usage: GWAStoolkit or2beta [options]\n\n"
    "Required:\n"
    "  --gwas-summary FILE\n"
    "  --out FILE\n"
    "  --or COL        OR column name\n"
    "  --SNP COL       SNP identifier column\n\n"
    "Required allele columns:\n"
    "  --A1 COL   (default: A1)\n"
    "  --A2 COL   (default: A2)\n\n"
    "Optional columns:\n"
    "  --se   COL   SE of OR (optional)\n"
    "  --pval COL   P-value (used to infer SE if OR_SE not provided)\n"
    "  --freq COL   (used in QC)\n"
    "  --n    COL   (used in QC)\n\n"
    "Optional:\n"
    "  --format gwas|cojo|popcorn|mrmega\n"
    "  --maf VAL\n"
    "  --remove-dup-snp\n"
    "  --threads N\n"
    "  --log FILE\n";
}

static void check_format_required_columns(const CommonArgs& P, const std::string& who)
{
    // gwas 格式不做严格限制，可仅输出原始列
    if (P.format == "gwas") return;
    // rsidImpu 不需要 SNP 列
    bool need_snp = (who != "rsidImpu");
    if (need_snp) {
        require(!P.col_SNP.empty(), who + " requires --SNP column for formatted output.");
    }

    require(!P.g_A1.empty(),    who + " requires --A1 column for formatted output.");
    require(!P.g_A2.empty(),    who + " requires --A2 column for formatted output.");
    require(!P.g_p.empty(),     who + " requires --pval column for formatted output.");
    require(!P.col_freq.empty(),who + " requires --freq column for formatted output.");
    require(!P.col_beta.empty(),who + " requires --beta column for formatted output.");
    require(!P.col_se.empty(),  who + " requires --se column for formatted output.");
    require(!P.col_n.empty(),   who + " requires --n column for formatted output.");
}

// ------------------------- 解析 rsid-impu -----------------------
Args_RsidImpu parse_args_rsidimpu(int argc, char* argv[]) {
    map<string,string> args;
    set<string> flags = {"--remove-dup-snp"};

    for (int i=1; i<argc; ) {
        string key = argv[i];

        if (key == "--help") {
            print_rsidimpu_help();
            exit(0);
        }

        // ★ unknown parameter check
        if (!common_params.count(key) &&
            !rsidimpu_params.count(key)) {
            LOG_ERROR("Unknown parameter: " + key);
            exit(1);
        }

        // flags
        if (flags.count(key)) {
            args[key] = "1"; i++; continue;
        }

        if (i+1 >= argc) {
            LOG_ERROR("Missing value for " + key);
            exit(1);
        }

        args[key] = argv[i+1];
        i += 2;
    }

    Args_RsidImpu P;
    parse_common(P, args);

    // for gwas columns
    require(!P.g_chr.empty(), "rsidImpu requires --chr column.");
    require(!P.g_pos.empty(), "rsidImpu requires --pos column.");
    require(!P.g_A1.empty(),  "rsidImpu requires --A1 column.");
    require(!P.g_A2.empty(),  "rsidImpu requires --A2 column.");

    // Required for rsid-impu
    require(args.count("--dbsnp"), "Missing required: --dbsnp");
    require(args.count("--dbchr"), "Missing required: --dbchr");
    require(args.count("--dbpos"), "Missing required: --dbpos");
    require(args.count("--dbA1"),  "Missing required: --dbA1");
    require(args.count("--dbA2"),  "Missing required: --dbA2");
    require(args.count("--dbrsid"),"Missing required: --dbrsid");

    P.dbsnp_file = args["--dbsnp"];
    P.d_chr  = args["--dbchr"];
    P.d_pos  = args["--dbpos"];
    P.d_A1   = args["--dbA1"];
    P.d_A2   = args["--dbA2"];
    P.d_rsid = args["--dbrsid"];

    check_format_required_columns(P, "rsidImpu");
    return P;
}

// ------------------------- 解析 convert ------------------------------
Args_Convert parse_args_convert(int argc, char* argv[]) {
    map<string,string> args;

    for (int i=1; i<argc; ) {
        string key = argv[i];

        if (key == "--help") {
            print_convert_help();
            exit(0);
        }

        // ★ unknown parameter check
        if (!common_params.count(key) &&
            !convert_params.count(key)) {
            LOG_ERROR("Unknown parameter: " + key);
            exit(1);
        }

        if (i+1 >= argc) {
            LOG_ERROR("Missing value for " + key);
            exit(1);
        }

        args[key] = argv[i+1];
        i += 2;
    }

    Args_Convert P;
    parse_common(P, args);

    check_format_required_columns(P, "convert");
    return P;
}

// ------------------------- 解析 or2beta ------------------------------
Args_Or2Beta parse_args_or2beta(int argc, char* argv[]) {
    map<string,string> args;

    for (int i=1; i<argc; ) {
        string key = argv[i];

        if (key == "--help") {
            print_or2beta_help();
            exit(0);
        }

        // ★ unknown parameter check
        if (!common_params.count(key) &&
            !or2beta_params.count(key)) {
            LOG_ERROR("Unknown parameter: " + key);
            exit(1);
        }

        if (i+1 >= argc) {
            LOG_ERROR("Missing value for " + key);
            exit(1);
        }

        args[key] = argv[i+1];
        i += 2;
    }

    Args_Or2Beta P;
    parse_common(P, args);

    require(args.count("--or"), "Missing required: --or");
    P.col_or = args["--or"];
    require(!P.col_SNP.empty(), "or2beta requires --SNP column.");
    require(!P.col_or.empty(),  "or2beta requires --or column.");

    check_format_required_columns(P, "or2beta");
    return P;
}