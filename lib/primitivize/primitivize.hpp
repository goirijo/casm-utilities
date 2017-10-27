#ifndef CASM_UTILS_PRIMITIVIZE
#define CASM_UTILS_PRIMITIVIZE

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string primitivize_launcher_name();

    void primitivize_utility_launch(int argc, char *argv[]);

    LaunchRuleList primitivize_initializer(po::options_description &primitivize_desc);
}


#endif
