#include "utils/args.hpp"
#include "utils/log.hpp"
#include "utils/gadgets.hpp"
#include "or2beta/or2beta.hpp"

int cmd_or2beta(int argc, char* argv[]) {

    Args_Or2Beta P = parse_args_or2beta(argc, argv);

    Gadget::Timer timer;
    timer.setTime();

    LOG_INFO("Running or2beta ...");
    run_or2beta(P);
    LOG_INFO("or2beta finished.");

    return 0;
}
