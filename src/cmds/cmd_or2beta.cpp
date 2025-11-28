#include "utils/args.hpp"
#include "utils/log.hpp"
#include "or2beta/or2beta.hpp"

#include <fstream>

int cmd_or2beta(int argc, char* argv[]) {

    Args_Or2Beta P = parse_args_or2beta(argc, argv);

    // log
    std::ofstream log_ofs;
    if (P.log_enabled) {
        log_ofs.open(P.log_file);
        g_log = &log_ofs;
    }
    g_log_to_console = true;

    // LOG_INFO("Running or2beta ... ");

    run_or2beta(P);

    // LOG_INFO("or2beta finished.");
    return 0;
}
