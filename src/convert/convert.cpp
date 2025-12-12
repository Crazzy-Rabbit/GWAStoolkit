#include "convert/convert.hpp"

#include "utils/linereader.hpp"
#include "utils/writer.hpp"
#include "utils/util.hpp"
#include "utils/log.hpp"
#include "utils/gwasQC.hpp"
#include "utils/FormatEngine.hpp"

#include <unordered_map>
#include <string>
#include <deque>
#include <algorithm>

using namespace std;

void run_convert(const Args_Convert& P){
    LineReader lr(P.gwas_file);
    string line;

    // read header
    if (!lr.getline(line)){
        LOG_ERROR("Empty GWAS summary file in convert.");
        exit(1);
    }
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    auto header = split(line);
    
    // header check
    int idx_snp  = find_col(header, P.col_SNP);
    require(idx_snp >= 0, "GWAS missing required column [" + P.col_SNP + "] for convert.");

    int idx_A1   = find_col(header, P.g_A1);
    require(idx_A1 >= 0, "GWAS missing required column [" + P.g_A1 + "] for convert.");

    int idx_A2   = find_col(header, P.g_A2);
    require(idx_A2 >= 0, "GWAS missing required column [" + P.g_A2 + "] for convert.");

    int idx_freq = find_col(header, P.col_freq);
    require(idx_freq >= 0, "GWAS missing required column [" + P.col_freq + "] for convert.");

    int idx_beta = find_col(header, P.col_beta);
    require(idx_beta >= 0, "GWAS missing required column [" + P.col_beta + "] for convert.");

    int idx_se   = find_col(header, P.col_se);
    require(idx_se >= 0, "GWAS missing required column [" + P.col_se + "] for convert.");

    int idx_p    = find_col(header, P.g_p);
    require(idx_p >= 0, "GWAS missing required column [" + P.g_p + "] for convert.");

    int idx_n    = find_col(header, P.col_n);
    require(idx_n >= 0, "GWAS missing required column [" + P.col_n + "] for convert.");

    // read all line
    deque<string> lines;
    while (lr.getline(line)) {
        if (!line.empty())
            lines.push_back(line);
    }
    size_t n = lines.size();
    LOG_INFO("Loaded GWAS lines for convert: " + to_string(n));

    // QC, default maf = 0.01
    vector<bool> keep(n,true);
    // 如果列都存在，就执行QC，否则给warning
    bool can_qc = (idx_beta>=0 || idx_se>=0 || idx_freq>=0 || idx_p>=0 || idx_n>=0);
    if (can_qc) {
        LOG_INFO("QC applied in partial-column mode.");
        gwas_basic_qc(lines, header,
                    idx_beta, idx_se, idx_freq, idx_p, idx_n,
                    keep, P.maf_threshold);
    } else {
        LOG_WARN("Cannot perform full QC in convert (missing beta/se/freq/N/P columns).");
    }

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

    // out format
    FormatEngine FE;
    FormatSpec spec = FE.get_format(P.format);
    Writer fout(P.out_file, P.format);

    if (!fout.good()){
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

    // write context
    for (size_t i=0; i<n; i++){
        if (!keep[i]) continue;

        lines[i].erase(remove(lines[i].begin(), lines[i].end(), '\r'), lines[i].end());
        auto f = split(lines[i]);

        // ★ 检查列数是否与 header 一致，否则跳过
        if ((int)f.size() != (int)header.size()) continue;
        if (f[idx_snp].empty()) continue;  

        if (P.format == "gwas"){
            fout.write_line(lines[i]);
            continue;
        }

        unordered_map<string,string> row;
        row["SNP"]  = f[idx_snp];
        row["A1"]   = f[idx_A1]; 
        row["A2"]   = f[idx_A2];
        row["freq"] = f[idx_freq];
        row["beta"] = f[idx_beta];
        row["se"]   = f[idx_se];
        row["p"]    = f[idx_p];
        row["N"]    = f[idx_n];

        string out_line = FE.format_line(spec, row);
        fout.write_line(out_line);
    }

    LOG_INFO("convert finished (format=" + P.format + ").");
}
