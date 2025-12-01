#include "computeNeff.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/FormatEngine.hpp"
#include "utils/gadgets.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

// Neff
static inline double calc_neff(double cs, double ct){
    double s = cs + ct;
    if (s <= 0) return 0;
    return 4.0 * cs * ct / s;
}

// standardization of beta abd se
static inline bool std_effect(
    double freq, double beta_old, double se_old, double Neff,
    double &beta_new, double &se_new
) {
    if (freq <= 0.0 || freq >= 1.0) return false;
    if (se_old <= 0.0) return false;
    if (!std::isfinite(Neff) || Neff <= 0.0) return false;

    double z = beta_old / se_old;
    double denom = 2.0 * freq * (1.0-freq) * (Neff + z*z);
    if (denom <= 0.0) return false;

    se_new   = 1.0 / sqrt(denom);
    beta_new = z * se_new;
    return true;
}

void run_computeNeff(const Args_CalNeff& P)
{

    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);

    LOG_INFO("computeNeff: loading GWAS...");
    Writer fout(P.out_file, P.format);
    if (!fout.good()){
        LOG_ERROR("Cannot open output: " + P.out_file);
        return;
    }

    LineReader reader(P.gwas_file);
    string line;
    if (!reader.getline(line)) {
        LOG_ERROR("Empty GWAS file: " + P.gwas_file);
        exit(1);
    }
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);

    // idx
    int idx_snp  = find_col(header, P.col_SNP);
    int idx_A1   = find_col(header, P.g_A1);
    int idx_A2   = find_col(header, P.g_A2);
    int idx_freq = find_col(header, P.col_freq);
    int idx_beta = find_col(header, P.col_beta);
    int idx_se   = find_col(header, P.col_se);
    int idx_p    = find_col(header, P.g_p);

    int idx_case    = -1;
    int idx_control = -1;

    // -------------------------
    // Case 1: per-SNP NEFF 计算
    // -------------------------
    if (P.is_column){
        idx_case    = find_col(header, P.case_col);
        idx_control = find_col(header, P.control_col);
        if (idx_case < 0) {
            LOG_ERROR("Cannot find case column: " + P.case_col);
            exit(1);
        }
        if (idx_control < 0) {
            LOG_ERROR("Cannot find control column: " + P.control_col);
            exit(1);
        }
    }

    // -------------------------
    // Case 2: fixed case/control
    // -------------------------
    double Neff_fixed = std::numeric_limits<double>::quiet_NaN();
    if (P.is_single){
        Neff_fixed = calc_neff(P.case_n, P.control_n);
        if (!std::isfinite(Neff_fixed) || Neff_fixed <= 0.0){
            LOG_ERROR("Invalid fixed case/control: " + 
                std::to_string(P.case_n) + "," + std::to_string(P.control_n));
            exit(1);
        }
        LOG_INFO("Fixed-mode Neff = " + std::to_string(Neff_fixed));

    }

    // write header
    {
        if (!spec.cols.empty()){
            string header_line;
            for (size_t i=0; i<spec.cols.size(); i++){
                if (i) header_line += "\t";
                header_line += spec.cols[i];
            }
            fout.write_line(header_line);
        }
    }

    // process GWAS per row
    while (reader.getline(line)){
        if (line.empty()) continue;
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        auto f = split(line);

        // 决定当前 SNP 的 Neff（根据两种模式）
        double Neff = NAN;
        if (P.is_single){
            Neff = Neff_fixed;
        } else if (P.is_column){
            if ((int)f.size() <= std::max(idx_case, idx_control)) continue;
            double cs = atof(f[idx_case].c_str());
            double ct = atof(f[idx_control].c_str());
            Neff      = calc_neff(cs, ct);
        }

        if (!std::isfinite(Neff) || Neff <= 0.0) {
            continue;
        }

        unordered_map<string,string> row;

        row["SNP"]  = f[idx_snp]; 
        row["A1"]   = f[idx_A1];
        row["A2"]   = f[idx_A2];

        // freq/beta/se scaling
        if (idx_freq>=0 && idx_beta>=0 && idx_se>=0 &&
            idx_freq < (int)f.size() &&
            idx_beta < (int)f.size() &&
            idx_se   < (int)f.size())
        {
            double freq_old = atof(f[idx_freq].c_str());
            double beta_old = atof(f[idx_beta].c_str());
            double se_old   = atof(f[idx_se].c_str());

            double beta_new, se_new;
            if (std_effect(freq_old, beta_old, se_old, Neff, beta_new, se_new)){
                row["freq"] = f[idx_freq];
                row["beta"] = std::to_string(beta_new);
                row["se"]   = std::to_string(se_new);
            } else {
                row["freq"] = f[idx_freq];
                row["beta"] = f[idx_beta];
                row["se"]   = f[idx_se];
            }
        }

        row["p"]    = f[idx_p];
        row["N"]    = std::to_string(Neff);
        fout.write_line(FE.format_line(spec, row));
    }
}