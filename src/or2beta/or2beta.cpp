#include "or2beta/or2beta.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/gwasQC.hpp"
#include "utils/FormatEngine.hpp"

#include <algorithm>
#include <unordered_map>
#include <cmath>

using namespace std;

// --------- P 值 → Z 值（gsMap 方法）------------
static double p_to_z(double p);

#ifdef USE_RMATH
    #include <Rmath.h>
    static double p_to_z(double p){
        if (p <= 0) return 38.0;
        if (p >= 1) return 0.0;
        return fabs(qnorm(p/2.0, 0.0, 1.0, 1, 0));
    }
#endif

void run_or2beta(const Args_Or2Beta& P){
    LOG_INFO("Running or2beta ...");

    LineReader lr(P.gwas_file);
    string line;

    if (!lr.getline(line)) {
        LOG_ERROR("Empty GWAS summary file in or2beta.");
        exit(1);
    }
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);

    // ---- locate columns (support global names) ----
    int idx_snp = find_col(header, P.col_SNP);
    int idx_A1  = find_col(header, P.g_A1);
    int idx_A2  = find_col(header, P.g_A2);
    int idx_or  = find_col(header, P.col_or);
    int idx_se  = find_col(header, P.col_se);
    int idx_p   = find_col(header, P.g_p);   // optional; if not found, p=1.0

    if (idx_snp<0 || idx_A1<0 || idx_A2<0 || idx_or<0 || idx_p<0) {
        LOG_ERROR("Missing required columns for or2beta: SNP/A1/A2/OR/P");
        exit(1);
    }

    // ---------- read lines ----------
    deque<string> lines;
    while (lr.getline(line)){
        if (!line.empty())
            lines.push_back(line);
    }

    size_t n = lines.size();
    LOG_INFO("Loaded " + to_string(n) + " GWAS lines for or2beta.");


    // ======================= QC（必须先做） =======================
    vector<bool> keep_qc(n, true);

    int idx_freq = find_col(header, P.col_freq);
    int idx_beta = -1;
    int idx_n    = find_col(header, P.col_n);

    bool can_qc = (idx_freq>=0 && idx_p>=0 && idx_n>=0);

    if (can_qc) {
        gwas_basic_qc(lines, header, 
            idx_beta, idx_se, idx_freq, idx_p, idx_n,
            keep_qc,
            P.maf_threshold);
    } else {
        LOG_WARN("QC not fully applied: missing freq/p/N columns.");
    }

    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);

    Writer fout(P.out_file, P.format);
    if (!fout.good()) {
        LOG_ERROR("Cannot open output file: " + P.out_file);
        exit(1);
    }

    // 写 header
    {
        string h;
        for (size_t i=0; i<spec.cols.size(); i++){
            if (i) h += "\t";
            h += spec.cols[i];
        }
        fout.write_line(h);
    }

    // process lines
    for (auto &ln : lines) {

        string line = ln;
        if (line.empty()) continue;
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        auto f = split(line);

        if ((int)f.size() <= idx_or) continue;

        double OR  = stod(f[idx_or]);
        double beta = log(OR);
        
        // ---- compute se
        double se   = NAN;
        if (idx_p >= 0){
            double pval = stod(f[idx_p]);
            double z = p_to_z(pval);
            se = (z > 0 ? fabs(beta) / z : 999);
        } else {
            se = 999;
        }

        // construct row for FormatEngine
        unordered_map<string,string> row;

        for (auto &col : spec.cols) {

            if (col == "SNP"   && idx_snp>=0) row[col] = f[idx_snp];
            else if (col == "A1" && idx_A1>=0) row[col] = f[idx_A1];
            else if (col == "A2" && idx_A2>=0) row[col] = f[idx_A2];
            else if (col == "beta") row[col] = to_string(beta);
            else if (col == "se")   row[col] = to_string(se);
            else if (col == "p" && idx_p>=0) row[col] = f[idx_p];
            else if (col == "P" && idx_p>=0) row[col] = f[idx_p];
            else
                row[col] = "";   // not available
        }

        fout.write_line(FE.format_line(spec, row));
    }

    LOG_INFO("or2beta finished.");
}