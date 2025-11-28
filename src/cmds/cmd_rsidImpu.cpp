#include "rsidImpu/dbsnp.hpp"
#include "rsidImpu/gwas.hpp"

#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"

#include <fstream>

int cmd_rsidImpu(int argc, char* argv[]){
    Args_RsidImpu P = parse_args_rsidimpu(argc, argv);

    // 初始化日志
    std::ofstream log_ofs;
    if (P.log_enabled) {
        log_ofs.open(P.log_file);
        g_log = &log_ofs;
    }

    g_log_to_console = true;

    Gadget::Timer timer;
    timer.setTime();

    LOG_INFO("Running rsidImpu ... ");

    auto db = load_dbsnp(P);
    process_gwas(P, db);

    LOG_INFO("rsidImpu finished.");
    LOG_INFO("Total time: " + timer.format(timer.getElapse()));

    return 0;
}