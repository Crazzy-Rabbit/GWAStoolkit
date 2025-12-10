#include "computeNeff.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/FormatEngine.hpp"
#include "utils/gadgets.hpp"
#include "utils/gwasQC.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <deque>

using namespace std;

// Neff
static inline double calc_neff(double cs, double ct){
    double s = cs + ct;
    if (s <= 0) return 0;
    return 4.0 * cs * ct / s;
}

// standardization of beta and se
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
    Writer fout(P.out_file, P.format);

    if (!fout.good()){
        LOG_ERROR("Cannot open output: " + P.out_file);
        return;
    }

    // header
    LineReader reader(P.gwas_file);
    string line;
    if (!reader.getline(line)) {
        LOG_ERROR("Empty GWAS file: " + P.gwas_file);
        exit(1);
    }
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);

    bool has_N = false;
    int idx_N  = -1;

    for (size_t i = 0; i < header.size(); i++){
        string hlow = header[i];
        std::transform(hlow.begin(), hlow.end(), hlow.begin(), ::tolower);

        if (hlow == "n"){
            has_N = true;
            idx_N = i;
            break;
        }
    }

    // header check
    int idx_snp = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for computeNeff.");

    int idx_A1 = find_col(header, P.g_A1);
    require(idx_A1 >= 0, "GWAS missing required column [" + P.g_A1 + "] for computeNeff.");

    int idx_A2 = find_col(header, P.g_A2);
    require(idx_A2 >= 0, "GWAS missing required column [" + P.g_A2 + "] for computeNeff.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for computeNeff.");

    int idx_beta = find_col(header, P.col_beta);
    require(idx_beta >= 0, "GWAS missing required column [" + P.col_beta + "] for computeNeff.");

    int idx_se = find_col(header, P.col_se);
    require(idx_se >= 0, "GWAS missing required column [" + P.col_se + "] for computeNeff.");

    int idx_p = find_col(header, P.g_p);
    require(idx_p >= 0, "GWAS missing required column [" + P.g_p + "] for computeNeff.");

    int idx_case    = -1;
    int idx_control = -1;

    // 读入所有数据行
    deque<string> lines;
    while (reader.getline(line)){
        if (!line.empty()) lines.push_back(line);
    }
    size_t n = lines.size();
    LOG_INFO("Loaded " + to_string(n) + " GWAS lines for computeNeff.");

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

    //================ QC + 去重 (按 SNP) ================
    vector<bool> keep(n, true);
    int idx_n = -1;
    
    bool can_qc = (idx_beta>=0 && idx_se>=0 && idx_freq>=0 && idx_p>=0);
    if (can_qc){
        int idx_n_safe = (idx_n >= 0 ? idx_n : idx_freq); 

        gwas_basic_qc(lines, header,
                    idx_beta, idx_se, idx_freq, idx_p, idx_n_safe,
                    keep, P.maf_threshold);
    } else {
        LOG_WARN("Cannot perform full QC in computeNeff (missing beta/se/freq/N/p columns).");
    }

    // remove dup SNP（用 SNP 作为 key）
    if (P.remove_dup_snp){
        vector<string> snp_vec(n);
        for (size_t i=0; i<n; i++){
            if (!keep[i]) continue;
            auto f = split(lines[i]);
            if (idx_snp>=0 && idx_snp<(int)f.size())
                snp_vec[i] = f[idx_snp];
        }
        gwas_remove_dup(lines, header, idx_p, snp_vec, keep);
    }


    // writer header
    if (P.format == "gwas") {
        // raw header
        string h;
        for (size_t i=0; i<header.size(); i++){
            if (i) h += "\t";
            h += header[i];
        }

        if (!has_N){
            h += "\tN";
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

    // process GWAS per row
    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;

        string ln = lines[i];
        ln.erase(remove(ln.begin(), ln.end(), '\r'), ln.end());
        auto f = split(ln);

        // ★ 检查列数是否与 header 一致，否则跳过
        if ((int)f.size() != (int)header.size()) continue;
        if (f[idx_snp].empty()) continue;  

        // 计算当前 SNP 的 Neff
        double Neff = NAN;
        if (P.is_single){
            Neff = Neff_fixed;
        } else if (P.is_column){
            if ((int)f.size() <= std::max(idx_case, idx_control)) continue;
            double cs = atof(f[idx_case].c_str());
            double ct = atof(f[idx_control].c_str());
            Neff      = calc_neff(cs, ct);
        }

        if (!std::isfinite(Neff) || Neff <= 0.0) continue;

        if (P.format == "gwas"){
            auto f = split(ln);

            if (has_N){
                // ★ 覆盖原 N
                if (idx_N < (int)f.size()){
                    f[idx_N] = std::to_string(Neff);
                }
                // 重建一行
                string out;
                for (size_t j=0; j<f.size(); j++){
                    if (j) out +="\t";
                    out += f[j];
                }
                fout.write_line(out);
            } else {
                // ★ 原来没有 N → 追加
                string out = ln + "\t" + std::to_string(Neff);
                fout.write_line(out);
            }
            continue;
        }

        unordered_map<string,string> row;

        row["SNP"] = f[idx_snp]; 
        row["A1"]  = f[idx_A1];
        row["A2"]  = f[idx_A2];

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

        if (idx_p >= 0 && idx_p < (int)f.size()) row["p"] = f[idx_p];
        row["N"] = std::to_string(Neff);

        string out_line = FE.format_line(spec, row);
        fout.write_line(out_line);
    }

}