//
//  args.cpp
//  GWAStookit
//
//  Created by Lulu Shi on 02/12/2025.
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
    "--dbsnp", "--dbchr", "--dbpos", "--dbA1", "--dbA2", "--dbrsid",
    "--chr", "--pos"
};
static const std::set<std::string> convert_params = {};
static const std::set<std::string> or2beta_params = {
    "--or"
};
static const set<string> calneff_params = {
    "--case", "--control",       // fixed-mode
    "--case-col", "--control-col" // per-SNP mode
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
    require(args.count("--out"),          "Missing required: --out");

    C.gwas_file = args["--gwas-summary"];
    C.out_file  = args["--out"];

    if (args.count("--threads"))
        C.threads = stoi(args["--threads"]);

    if (args.count("--log")) {
        C.log_enabled = true;
        C.log_file    = args["--log"];
    }

    if (args.count("--remove-dup-snp")){
        C.remove_dup_snp = true;
    }
    
    if (args.count("--maf")) {
        C.maf_threshold = stod(args["--maf"]);
    }

    // 公共列名
    if (args.count("--SNP"))   C.col_SNP  = args["--SNP"];
    if (C.col_SNP.empty())     C.col_SNP  = "SNP";

    if (args.count("--chr"))   C.g_chr    = args["--chr"];
    if (C.g_chr.empty())       C.g_chr    = "CHR";

    if (args.count("--pos"))   C.g_pos    = args["--pos"];
    if (C.g_pos.empty())       C.g_pos    = "POS";

    if (args.count("--A1"))    C.g_A1     = args["--A1"];
    if (C.g_A1.empty())        C.g_A1     = "A1";

    if (args.count("--A2"))    C.g_A2     = args["--A2"];
    if (C.g_A2.empty())        C.g_A2     = "A2";

    if (args.count("--pval"))  C.g_p      = args["--pval"];
    if (C.g_p.empty())         C.g_p      = "p";

    if (args.count("--freq"))  C.col_freq = args["--freq"];
    if (C.col_freq.empty())    C.col_freq = "freq";

    if (args.count("--beta"))  C.col_beta = args["--beta"];
    if (C.col_beta.empty())    C.col_beta = "b";

    if (args.count("--se"))    C.col_se   = args["--se"];
    if (C.col_se.empty())      C.col_se   = "se";

    if (args.count("--n"))     C.col_n    = args["--n"];
    if (C.col_n.empty())       C.col_n    = "N";

    if (args.count("--format")) 
        C.format = args["--format"];
    else
        C.format = "gwas";
    check_format(C.format);
}

// ======================================================
//                     HELP 信息
// ======================================================
void print_rsidimpu_help() {
    cerr <<
    "Usage:\n"
    "  GWAStoolkit rsidImpu [options]\n\n"

    "Description:\n"
    "  Annotate GWAS summary statistics with dbSNP rsID.\n"
    "  Allele matching supports flips and strand complements.\n\n"

    "Required arguments:\n"
    "  --gwas-summary FILE        Input GWAS summary statistics (txt / tsv / gz)\n"
    "  --dbsnp FILE               dbSNP or PLINK .bim file (txt / gz)\n"
    "  --out FILE                 Output file (txt or .gz)\n"

    "Required dbSNP columns:\n"
    "  --dbchr  COL  Chromosome column      (default: CHR)\n"
    "  --dbpos  COL  Base position column   (default: POS)\n"
    "  --dbrsid COL  rsid for SNP           (default: ID)\n"
    "  --dbA1   COL  REF allele             (default: REF)\n"
    "  --dbA2   COL  ALT allele             (default: ALT)\n\n" 

    "Required GWAS columns:\n"
    "  --chr  COL   Chromosome column      (default: CHR)\n"
    "  --pos  COL   Base position column   (default: POS)\n"
    "  --A1   COL   Effect allele          (default: A1)\n"
    "  --A2   COL   Other allele           (default: A2)\n\n"

    "Optional output format:\n"
    "  --format gwas|cojo|popcorn|mrmega   (default: gwas)\n\n"

    "Optional GWAS columns (required depending on --format):\n"
    "  --freq COL   Allele frequency       (default: freq)\n"
    "  --beta COL   Effect size            (default: b)\n"
    "  --se   COL   Standard error         (default: se)\n"
    "  --pval COL   P-value column         (default: p)\n"
    "  --n    COL   Sample size            (default: N)\n\n"

    "Quality Control options:\n"
    "  --maf VAL            MAF threshold (default: 0.01)\n"
    "  --remove-dup-snp     Keep only lowest-P SNP if duplicates exist\n\n"

    "Other options:\n"
    "  --threads N          Number of threads (default: 1)\n"
    "  --log FILE           Write log output to FILE\n";
}

void print_convert_help() {
    cerr <<
    "Usage:\n"
    "  GWAStoolkit convert [options]\n\n"

    "Description:\n"
    "  Convert GWAS summary statistics into specific downstream formats.\n"
    "  Supported: gwas, cojo, popcorn, mrmega.\n\n"

    "Required arguments:\n"
    "  --gwas-summary FILE     Input GWAS summary statistics (txt / gz)\n"
    "  --out FILE              Output file (txt or .gz)\n"
    "  --format gwas|cojo|popcorn|mrmega\n"
    "  --SNP COL               SNP identifier column\n\n"

    "Required GWAS columns for conversion:\n"
    "  --SNP  COL   Marker name          (default: SNP)\n"
    "  --A1   COL   Effect allele        (default: A1)\n"
    "  --A2   COL   Other allele         (default: A2)\n"
    "  --freq COL   Allele frequency     (default: freq)\n"
    "  --beta COL   Beta                 (default: b)\n"
    "  --se   COL   Standard error       (default: se)\n"
    "  --pval COL   P-value              (default: p)\n"
    "  --n    COL   Sample size          (default: N)\n\n"

    "Quality Control options:\n"
    "  --maf VAL            MAF threshold (default: 0.01)\n"
    "  --remove-dup-snp     Keep only lowest-P SNP if duplicates exist\n\n"

    "Other options:\n"
    "  --threads N\n"
    "  --log FILE\n";
}

void print_or2beta_help() {
    cerr <<
    "Usage:\n"
    "  GWAStoolkit or2beta [options]\n\n"

    "Description:\n"
    "  Convert Odds Ratio (OR) to Beta and SE.\n"
    "  If SE missing, SE is inferred from p-value.\n\n"

    "Required arguments:\n"
    "  --gwas-summary FILE    Input GWAS summary statistics (txt / gz)\n"
    "  --out FILE             Output file (txt or .gz)\n"

    "Required GWAS columns for or2beta:\n"
    "  --SNP  COL   Marker name          (default: SNP)\n"
    "  --A1   COL   Effect allele        (default: A1)\n"
    "  --A2   COL   Other allele         (default: A2)\n"
    "  --freq COL   Allele frequency     (default: freq)\n"
    "  --or   COL   OR values            (default: OR)\n"
    "  --pval COL   P-value              (default: p)\n"

    "Quality Control options:\n"
    "  --maf VAL\n"
    "  --remove-dup-snp\n\n"

    "Optional output format:\n"
    "  --format gwas|cojo|popcorn|mrmega   (default: gwas)\n\n"

    "Other options:\n"
    "  --threads N\n"
    "  --log FILE\n";
}

void print_calneff_help() {
    cerr <<
    "Usage:\n"
    "  GWAStoolkit computeNeff [options]\n\n"

    "Description:\n"
    "  Compute effective sample size for binary trait:\n"
    "      Neff = 4 * case * control / (case + control)\n"
    "  Supports fixed case/control counts or per-SNP columns.\n\n"

    "Mode 1: Fixed case/control values:\n"
    "  --case INT --control INT --out FILE\n\n"

    "Mode 2: Per-SNP counts:\n"
    "  --gwas-summary FILE --case-col COL --control-col COL --out FILE\n\n"

    "Required GWAS columns for computeNeff:\n"
    "  --SNP  COL   Marker name          (default: SNP)\n"
    "  --A1   COL   Effect allele        (default: A1)\n"
    "  --A2   COL   Other allele         (default: A2)\n"
    "  --freq COL   Allele frequency     (default: freq)\n"
    "  --beta COL   Beta                 (default: b)\n"
    "  --se   COL   Standard error       (default: se)\n"
    "  --pval COL   P-value              (default: p)\n"

    "Quality Control options:\n"
    "  --maf VAL            MAF threshold (default: 0.01)\n"
    "  --remove-dup-snp     Keep only lowest-P SNP if duplicates exist\n\n"

    "Optional output format:\n"
    "  --format gwas|cojo|popcorn|mrmega   (default: gwas)\n\n"

    "Other options:\n"
    "  --threads N\n"
    "  --log FILE\n";
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

    // Required for rsid-impu
    require(args.count("--dbsnp"), "Missing required: --dbsnp");
    P.dbsnp_file = args["--dbsnp"];

    if (args.count("--dbchr")) P.d_chr = args["--dbchr"]; else P.d_chr = "CHR";
    if (args.count("--dbpos")) P.d_pos = args["--dbpos"]; else P.d_pos = "POS";
    if (args.count("--dbA1"))  P.d_A1  = args["--dbA1"];  else P.d_A1  = "REF";
    if (args.count("--dbA2"))  P.d_A2  = args["--dbA2"];  else P.d_A2  = "ALT";
    if (args.count("--dbrsid"))P.d_rsid= args["--dbrsid"];else P.d_rsid= "ID";

    return P;
}

// ------------------------- 解析 convert ------------------------------
Args_Convert parse_args_convert(int argc, char* argv[]) {
    map<string,string> args;
    set<string> flags = {"--remove-dup-snp"};

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

    Args_Convert P;
    parse_common(P, args);

    require(!P.col_SNP.empty(),  "convert requires --SNP column.");
    require(!P.g_A1.empty(),     "convert requires --A1 column.");
    require(!P.g_A2.empty(),     "convert requires --A2 column.");
    require(!P.col_freq.empty(), "convert requires --freq column.");
    require(!P.col_beta.empty(), "convert requires --beta column.");
    require(!P.col_se.empty(),   "convert requires --se column.");
    require(!P.g_p.empty(),      "convert requires --pval column.");
    require(!P.col_n.empty(),    "convert requires --n column.");

    return P;
}

// ------------------------- 解析 or2beta ------------------------------
Args_Or2Beta parse_args_or2beta(int argc, char* argv[]) {
    map<string,string> args;
    set<string> flags = {"--remove-dup-snp"};

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

    Args_Or2Beta P;
    parse_common(P, args);

    require(args.count("--or"), "Missing required: --or");
    P.col_or = args["--or"];

    if (args.count("--or")) P.col_or= args["--or"]; else P.col_or= "OR";
    require(!P.col_SNP.empty(),  "or2beta requires --SNP column.");
    require(!P.g_A1.empty(),     "or2beta requires --A1 column.");
    require(!P.g_A2.empty(),     "or2beta requires --A2 column.");
    require(!P.col_freq.empty(), "or2beta requires --freq column.");
    require(!P.col_se.empty() || !P.g_p.empty(),
            "or2beta requires --se or --pval to compute SE.");

    return P;
}

Args_CalNeff parse_args_calneff(int argc, char* argv[])
{
    map<string,string> args;
    set<string> flags = {"--remove-dup-snp"};

    for (int i=1; i<argc; ) {
        string key = argv[i];

        if (key == "--help") {
            print_calneff_help();
            exit(0);
        }

        if (!common_params.count(key) &&
            !calneff_params.count(key)) {
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

    Args_CalNeff P;
    parse_common(P, args);

    bool fixed = args.count("--case")      && args.count("--control");
    bool perSNP = args.count("--case-col") && args.count("--control-col");

    require(fixed || perSNP, 
        "computeNeff requires --case/--control OR --case-col/--control-col");
    require(!(fixed && perSNP), 
        "Cannot mix fixed and per-SNP modes.");

    if (fixed) {
        P.is_single = true;
        P.case_n    = stoi(args["--case"]);
        P.control_n = stoi(args["--control"]);
    }

    if (perSNP) {
        P.is_column   = true;
        P.case_col    = args["--case-col"];
        P.control_col = args["--control-col"];
    }

    require(!P.col_SNP.empty(),  "computeNeff requires --SNP column.");
    require(!P.col_freq.empty(), "computeNeff requires --freq column.");
    require(!P.col_beta.empty(), "computeNeff requires --beta column.");
    require(!P.col_se.empty(),   "computeNeff requires --se column.");
    require(!P.g_p.empty(),      "computeNeff requires --pval column.");
    require(!P.col_n.empty(),    "computeNeff requires --n column.");

    return P;
}