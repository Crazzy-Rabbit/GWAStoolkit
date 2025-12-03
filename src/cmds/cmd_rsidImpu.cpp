#include "rsidImpu/dbsnp.hpp"
#include "rsidImpu/rsidImpu.hpp"

#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"

int cmd_rsidImpu(int argc, char* argv[]){
    Args_RsidImpu P = parse_args_rsidimpu(argc, argv);

    Gadget::Timer timer;
    timer.setTime();

    LOG_INFO("Running rsidImpu ... ");
    process_rsidImpu(P);
    LOG_INFO("rsidImpu finished.");
    return 0;
}