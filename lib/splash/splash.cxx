#include <casm/CASM_global_definitions.hh>
#include <iostream>
#include <functional>

#include "lib/splash/splash.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"


namespace casmUtilities
{
    std::string splash_launcher_name()
    {
        return "splash";
    }

    LaunchRuleList splash_initializer(po::options_description &splash_desc)
    {
        utilityProgramOptions::add_help_suboption(splash_desc);

        splash_desc.add_options()
            ("print,p", "Print the CASM logo to the screen")
            ("dumb-print,d", "Print the CASM logo as a word")
            ("squelch,s", "Don't print anything")
            ("number,n", po::value<int>()->default_value(1), "How many times to print the logo");

        LaunchRuleList splash_rules;
        splash_rules.add_any_inclusion("number",std::vector<std::string>{"dumb-print","print"});
        //splash_rules.add_any_exclusion("squelch", std::vector<std::string>{"dumb-print","print"});

        return splash_rules;
    }

    void splash_utility_launch(int argc, char *argv[])
    {
        Launcher splash_launch(argc, argv, splash_launcher_name(), casmUtilities::splash_initializer);

        if(splash_launch.count("help"))
        {
            std::cout<<splash_launch.utility().desc()<<std::endl;
        }

        if(splash_launch.count("print"))
        {
            CASM::print_splash(std::cout);
        }

        if(splash_launch.count("dumb-print"))
        {
            std::cout<<"CASM"<<std::endl;
        }

        if(splash_launch.count("test"))
        {
            std::cout<<"You tested the value "<<splash_launch.fetch<int>("test")<<std::endl;
        }

        return;
    }
}
