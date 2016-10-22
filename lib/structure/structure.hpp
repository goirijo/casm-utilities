#ifndef CASM_UTILS_STRUCTURE
#define CASM_UTILS_STRUCTURE

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string structure_launcher_name();

    void structure_utility_launch(int argc, char *argv[]);

    LaunchRuleList structure_initializer(po::options_description &structure_desc);
}


#endif
