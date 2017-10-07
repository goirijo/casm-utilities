#ifndef CASM_UTILS_BOXIFY_HH
#define CASM_UTILS_BOXIFY_HH

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string boxify_launcher_name();

    LaunchRuleList boxify_initializer(po::options_description &boxify_desc);

    void boxify_utility_launch(int argc, char *argv[]);
}

#endif
