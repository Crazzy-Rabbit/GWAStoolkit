#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"
#include "convert/convert.hpp"

int cmd_convert(int argc, char* argv[]){
    Args_Convert P = parse_args_convert(argc, argv);

    Gadget::Timer timer;
    timer.setTime();

    LOG_INFO("Running convert ...");
    run_convert(P);
    // finish

    return 0;
}

