#include<iostream>
#include<string>

#include "lib/definitions.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"

namespace casmUtilities
{
    std::string tester_launcher_name()
    {
        return "tester";
    }
    
    LaunchRuleList tester_initializer(po::options_description &tester_desc)
    {
        utilityProgramOptions::add_help_suboption(tester_desc);

        tester_desc.add_options()
            ("atest,a", "Type A test")
            ("parrot,p", po::value<std::string>(), "Repeats what you say")
            ("creative-parrot,c", po::value<std::string>()->implicit_value("SKWAK!"), "Repeats what you say but has its own ideas")
            ("defaulted,d", po::value<std::string>()->default_value("ARR!"), "Always called");

        LaunchRuleList tester_rules;

        return tester_rules;
    }

    void tester_utility_launch(int argc, char *argv[])
    {
        Launcher tester_launch(argc, argv, tester_launcher_name(), tester_initializer);

        if(tester_launch.count("help"))
        {
            std::cout<<tester_launch.utility().desc()<<std::endl;
            return;
        }

        if(tester_launch.count("atest"))
        {
            std::cout<<"This is just A test, stay calm."<<std::endl;
        }

        if(tester_launch.count("parrot"))
        {
            std::cout<<tester_launch.fetch<std::string>("parrot")<<std::endl;
        }

        if(tester_launch.count("creative-parrot"))
        {
            std::cout<<tester_launch.fetch<std::string>("creative-parrot")<<std::endl;
        }

        if(tester_launch.count("defaulted"))
        {
            std::cout<<tester_launch.fetch<std::string>("defaulted")<<std::endl;
        }

        return;

    }
}
