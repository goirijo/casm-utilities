#include <casm/CASM_global_definitions.hh>
#include <iostream>

#include "lib/splash/splash.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"


namespace casmUtilities
{
    void splash_utility_launch(int argc, char *argv[])
    {
        Launcher<casmUtilities::SplashOption> splash_launch(argc, argv);

        //m_argument_rules.add_exclusion("squelch","dumb-print");
        //m_argument_rules.add_exclusion("squelch","print");

        if(splash_launch.count("help"))
        {
            std::cout<<splash_launch.option().desc()<<std::endl;
        }

        if(splash_launch.count("print"))
        {
            CASM::print_splash(std::cout);
        }

        if(splash_launch.count("test"))
        {
            std::cout<<"You tested the value "<<splash_launch.vm()["test"].as<int>()<<std::endl;
        }


        return;
    }
}
