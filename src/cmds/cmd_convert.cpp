#include "utils/args.hpp"
#include "utils/log.hpp"
#include "convert/convert.hpp"

#include <fstream>

int cmd_convet(int argc, char* argv[]){
    Args_Convert P = parse_args_convert(argc, argv);

    // 初始化日志（与 rsidImpu 一致）
    std::ofstream log_ofs;
    if (P.log_enabled) {
        log_ofs.open(P.log_file);
        if (!log_ofs) {
            LOG_ERROR("Cannot open log file: " + P.log_file);
            return 1;
        }
        g_log = &log_ofs;
    }
    g_log_to_console = true;

    // LOG_INFO("Running convert ...");

    run_convert(P);

    // LOG_INFO("convert finished.");
    return 0;
}