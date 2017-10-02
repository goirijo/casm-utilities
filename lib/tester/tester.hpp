#ifndef CASM_UTILS_TESTER_HH
#define CASM_UTILS_TESTER_HH

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string tester_launcher_name();

    LaunchRuleList tester_initializer(po::options_description &tester_desc);

    void tester_utility_launch(int argc, char *argv[]);
}

#endif
