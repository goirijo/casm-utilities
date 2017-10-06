#ifndef CASM_UTILS_BEDAZZLE_HH
#define CASM_UTILS_BEDAZZLE_HH

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string bedazzle_launcher_name();

    LaunchRuleList bedazzle_initializer(po::options_description &bedazzle_desc);

    void bedazzle_utility_launch(int argc, char *argv[]);
}

#endif
