#include "or2beta/or2beta.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/gwasQC.hpp"
#include "utils/FormatEngine.hpp"
#include "utils/StatFunc.hpp"

#include <algorithm>
#include <unordered_map>
#include <cmath>

using namespace std;

void run_or2beta(const Args_Or2Beta& P){
    LineReader lr(P.gwas_file);
    string line;

    if (!lr.getline(line)) {
        LOG_ERROR("Empty GWAS summary file in or2beta.");
        exit(1);
    }
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);

    // header check
    int idx_snp = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for or2beta.");

    int idx_A1  = find_col(header, P.g_A1);
    require(idx_A1 >= 0,  "GWAS missing required column [" + P.g_A1 + "] for or2beta.");

    int idx_A2  = find_col(header, P.g_A2);
    require(idx_A2 >= 0,  "GWAS missing required column [" + P.g_A2 + "] for or2beta.");

    int idx_or  = find_col(header, P.col_or);
    require(idx_or >= 0, "GWAS missing required column [" + P.col_or + "] for or2beta.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for or2beta.");

    int idx_se   = find_col(header, P.col_se);
    int idx_p    = find_col(header, P.g_p);
    require(idx_se >= 0 || idx_p >= 0,
            "or2beta requires either SE column [" + P.col_se +
            "] or P column [" + P.g_p + "].");

    int idx_n    = find_col(header, P.col_n);

    // ---------- read lines ----------
    deque<string> lines;
    while (lr.getline(line)){
        if (!line.empty()) lines.push_back(line);
    }

    size_t n = lines.size();
    LOG_INFO("Loaded " + to_string(n) + " GWAS lines for or2beta.");

    // ======================= QC =======================
    vector<bool> keep(n, true);

    bool can_qc = (idx_freq>=0 && idx_p>=0 && idx_n>=0);
    if (can_qc) {
        gwas_basic_qc(lines, header, 
            -1, idx_se, idx_freq, idx_p, idx_n,
            keep,
            P.maf_threshold);
    } else {
        LOG_WARN("QC not fully applied: missing freq/p/N columns.");
    }

    // remove dup SNP：用 SNP 作为 key
    if (P.remove_dup_snp) {
        vector<string> snp_vec(n);
        for (size_t i=0;i<n;i++){
            if (!keep[i]) continue;
            auto f = split(lines[i]);
            if (idx_snp>=0 && idx_snp<(int)f.size())
                snp_vec[i] = f[idx_snp];
        }
        gwas_remove_dup(lines, header, idx_p, snp_vec, keep);
    }

    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);
    Writer fout(P.out_file, P.format);

    if (!fout.good()) {
        LOG_ERROR("Cannot open output file: " + P.out_file);
        exit(1);
    }

    // writer header
    if (P.format == "gwas") {
        // raw header
        string h;
        for (size_t i=0; i<header.size(); i++){
            if (i) h += "\t";
            h += header[i];
        }
        fout.write_line(h);
    } else {
        string h;
        for (size_t i=0; i<spec.cols.size(); i++){
            if (i) h += "\t";
            h += spec.cols[i];
        }
        fout.write_line(h);
    }

    // process lines
    for (size_t i=0; i<n; i++) {
        if (!keep[i]) continue;

        string ln = lines[i];
        ln.erase(remove(ln.begin(), ln.end(), '\r'), ln.end());
        auto f = split(ln);

        // ★ 检查列数是否与 header 一致，否则跳过
        if ((int)f.size() != (int)header.size()) continue;
        if (f[idx_snp].empty()) continue;  

        if ((int)f.size() <= idx_or) continue;

        double OR  = stod(f[idx_or]);
        double beta = log(OR);
        
        // 计算 se
        double se = NAN;
        if (idx_se >= 0 && idx_se < (int)f.size()){
            se = atof(f[idx_se].c_str());
        } else if (idx_p >= 0 && idx_p < (int)f.size()){
            double pval = atof(f[idx_p].c_str());
            double z    = StatFunc::p2z_two_tailed(pval);
            se = (z > 0 ? fabs(beta) / z : 999);
        } else {
            se = 999;
        }

        unordered_map<string,string> row;
        row["SNP"]  = f[idx_snp];
        row["A1"]   = f[idx_A1];
        row["A2"]   = f[idx_A2];
        if (idx_freq>=0 && idx_freq < (int)f.size()) row["freq"] = f[idx_freq];
        if (idx_n>=0 && idx_n < (int)f.size())       row["N"]    = f[idx_n];
        if (idx_p>=0 && idx_p < (int)f.size())       row["p"]    = f[idx_p];
        row["beta"] = to_string(beta);
        row["se"]   = to_string(se);

        string out_line = FE.format_line(spec, row);
        fout.write_line(out_line);
    }

}