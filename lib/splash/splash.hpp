#ifndef CASM_UTILS_SPLASH
#define CASM_UTILS_SPLASH

#include "lib/launch/launchers.hpp"

namespace casmUtilities
{
    class LaunchRuleList;

    std::string splash_name();

    LaunchRuleList splash_initializer(po::options_description &splash_desc);

    void splash_utility_launch(int argc, char *argv[]);
}


#endif
