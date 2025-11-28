//

#include "utils/log.hpp"
#include <iostream>
#include <string>

using namespace std;

int cmd_rsidImpu(int argc, char* argv[]);
int cmd_convert(int argc, char* argv[]);
int cmd_or2beta(int argc, char* argv[]);

void print_main_help() {
    cout << "Available commands:\n"
        << "   rsidImpu       Annotate GWAS sumstats with rsid\n"
        << "   convert        Convert GWAS format (GWAS, COJO, SMR, LDSC, MR-MEGA)\n"
        << "   or2beta        Convert OR to beta and SE\n\n"
        << "Example:\n"
        << "  GWAStoolkit <command> [options]\n\n";
}

int main(int argc, char* argv[]) {
    cout << "******************************************************************\n";
    cout << "* GWAStoolkit                                                    *\n";
    cout << "* A tool to treat GWAS summary statistics                        *\n";
    cout << "* Authors: Loren Shi                                             *\n";
    cout << "* MIT License                                                    *\n";
    cout << "******************************************************************\n\n";
    
    if (argc < 2) {
        print_main_help();
        return 0;
    }

    string cmd = argv[1];

    if (cmd == "rsidImpu")
        return cmd_rsidImpu(argc - 1, argv + 1);
    if (cmd == "convert")
        return cmd_convert(argc - 1, argv + 1);
    if (cmd == "or2beta")
        return cmd_or2beta(argc - 1, argv + 1);
    
    LOG_ERROR("Unknown command: " + cmd);
    print_main_help();
    return 1;
}