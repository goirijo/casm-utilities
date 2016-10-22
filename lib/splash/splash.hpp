#ifndef CASM_UTILS_SPLASH
#define CASM_UTILS_SPLASH

#include "lib/definitions.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string splash_launcher_name();

    LaunchRuleList splash_initializer(po::options_description &splash_desc);

    void splash_utility_launch(int argc, char *argv[]);
}


#endif
