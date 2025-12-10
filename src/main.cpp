#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"

#include <iostream>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int cmd_rsidImpu(int argc, char* argv[]);
int cmd_convert(int argc, char* argv[]);
int cmd_or2beta(int argc, char* argv[]);
int cmd_computeNeff(int argc, char* argv[]);

void print_main_help() {
    cerr << "Available commands:\n"
        << "   rsidImpu       Annotate GWAS sumstats with rsid\n"
        << "   convert        Convert GWAS format (GWAS, COJO, SMR, LDSC, MR-MEGA)\n"
        << "   or2beta        Convert OR to beta and SE\n"
        << "   computeNeff    Compute effect sample size for binary traits\n\n"
        << "Example:\n"
        << "  GWAStoolkit <command> [options]\n\n";
}

int main(int argc, char* argv[]) {
    cout << "************************************************\n";
    cout << "* GWAStoolkit                                  *\n";
    cout << "* A tool to treat GWAS summary statistics      *\n";
    cout << "* Authors: Loren Shi                           *\n";
    cout << "* MIT License                                  *\n";
    cout << "************************************************\n\n";
    
    if (argc < 2 || 
        std::string(argv[1]) == "--help" ||
        std::string(argv[1]) == "-h") {
        print_main_help();
        return 0;
    }

    string cmd = argv[1];

    int threads = 1;
    {
        // 在所有命令的 args 解析前读取 threads 参数
        for (int i=2; i<argc; i++){
            std::string x = argv[i];
            if (x == "--threads" && i+1 < argc) {
                threads = std::stoi(argv[i+1]);
            }
        }
    }

#ifdef _OPENMP
    if (threads > 0){
        omp_set_num_threads(threads);
        std::cout << "[INFO] Using threads = " << threads << "\n";
    }
#endif

    for (int i=2; i<argc; i++){
        if (std::string(argv[i]) == "--log" && i+1 < argc) {
            static std::ofstream log_ofs(argv[i+1]);
            if (!log_ofs){
                LOG_ERROR("ERROR: cannot open log file");
                return 1;
            }
            g_log = &log_ofs;
        }
    }
    g_log_to_console = true;

    // Timer
    Gadget::Timer timer;
    timer.setTime();
    
    LOG_INFO(string("Analysis started: ") + timer.getDate());

    int ret = 0;
    if (cmd == "rsidImpu") {
        ret = cmd_rsidImpu(argc-1, argv+1);
    }
    else if (cmd == "convert") {
        ret = cmd_convert(argc-1, argv+1);
    }
    else if (cmd == "or2beta") {
        ret = cmd_or2beta(argc-1, argv+1);
    }
    else if (cmd == "computeNeff") {
        ret = cmd_computeNeff(argc-1, argv+1);
    }
    else {
        LOG_ERROR("Unknown command: " + cmd);
        return 1;
    }

    timer.getTime();
    LOG_INFO(string("Analysis finished: ") + timer.getDate());
    LOG_INFO(string("Total runtime: ") + timer.format(timer.getElapse()));
    return ret;
}