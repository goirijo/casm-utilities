#include <casm/completer/Complete.hh>
#include <casm/completer/Handlers.hh>
#include <iostream>

#include "lib/completer/handlers.hpp"
#include "lib/splash/splash.hpp"
#include "lib/structure/structure.hpp"
#include "lib/tester/tester.hpp"
#include "lib/bedazzle/bedazzle.hpp"

namespace casmUtilities
{
    CASM::Completer::Engine build_casm_utilities_engine()
    {
        CASM::Completer::Engine casm_utilities_engine;

        UtilityHandler dumbsplash(splash_launcher_name(), casmUtilities::splash_initializer);
        casm_utilities_engine.push_back(CASM::Completer::Option(dumbsplash.tag(), dumbsplash.desc()));

        UtilityHandler dumbstructure(structure_launcher_name(), casmUtilities::structure_initializer);
        casm_utilities_engine.push_back(CASM::Completer::Option(dumbstructure.tag(), dumbstructure.desc()));

        UtilityHandler dumbtester(tester_launcher_name(), casmUtilities::tester_initializer);
        casm_utilities_engine.push_back(CASM::Completer::Option(dumbtester.tag(), dumbtester.desc()));

        UtilityHandler dumbbedazzle(bedazzle_launcher_name(), casmUtilities::bedazzle_initializer);
        casm_utilities_engine.push_back(CASM::Completer::Option(dumbbedazzle.tag(), dumbbedazzle.desc()));

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
