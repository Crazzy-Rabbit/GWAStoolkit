#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"
#include "computeNeff/computeNeff.hpp"

int cmd_computeNeff(int argc, char* argv[])
{
    Args_CalNeff P = parse_args_calneff(argc, argv);

    Gadget::Timer timer;
    timer.setTime();

    LOG_INFO("Running computeNeff ...");
    run_computeNeff(P);
    LOG_INFO("computeNeff finished.");

    return 0;
}