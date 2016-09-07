#include <casm/completer/Complete.hh>
#include <casm/completer/Handlers.hh>
#include <iostream>

#include "lib/completer/handlers.hpp"

namespace casmUtilities
{
    CASM::Completer::Engine build_casm_utilities_engine()
    {
        CASM::Completer::Engine casm_utilities_engine;

        SplashOption dumbsplash;
        casm_utilities_engine.push_back(CASM::Completer::Option(dumbsplash.tag(), dumbsplash.desc()));

        return casm_utilities_engine;
    }
}


using namespace casmUtilities;

int main(int argc, char *argv[])
{
    CASM::Completer::Engine casm_utilities_engine=build_casm_utilities_engine();


    if(argc == 1)
    {
        std::cout << casm_utilities_engine.probe_options();
    }

    if(argc == 2)
    {
        std::cout << casm_utilities_engine.probe_suboptions(argv[1]);
    }

    if(argc == 3)
    {
        std::cout << casm_utilities_engine.probe_arguments(argv[1], argv[2]);
    }


    return 0;
}
