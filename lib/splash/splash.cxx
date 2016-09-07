#include <iostream>
#include <functional>

#include "lib/splash/splash.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"


namespace casmUtilities
{
    std::string splash_name()
    {
        return "splash";
    }

    LaunchRuleList splash_initializer(po::options_description &splash_desc)
    {
        LaunchRuleList splash_rules;

        splash_desc.add_options()
            ("print,p", "Print the CASM logo to the screen")
            ("dumb-print,d", "Print the CASM logo as a word")
            ("squelch,s", "Don't print anything")
            ("test,t", CASM::po::value<int>(), "Test without storage")
            ("number,n", CASM::po::value<int>()->default_value(1), "How many times to print the logo");

        return splash_rules;
    }

    void splash_utility_launch(int argc, char *argv[])
    {
        Launcher splash_launch(argc, argv, splash_name(), casmUtilities::splash_initializer);

        if(splash_launch.count("help"))
        {
            std::cout<<splash_launch.utility().desc()<<std::endl;
        }

        if(splash_launch.count("print"))
        {
            CASM::print_splash(std::cout);
        }

        if(splash_launch.count("test"))
        {
            std::cout<<"You tested the value "<<splash_launch.fetch<int>("test")<<std::endl;
        }

        return;
    }
}
