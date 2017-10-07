#include <iostream>

#include "lib/completer/handlers.hpp"
#include "lib/splash/splash.hpp"
#include "lib/structure/structure.hpp"
#include "lib/tester/tester.hpp"
#include "lib/bedazzle/bedazzle.hpp"
#include "lib/boxify/boxify.hpp"

using namespace casmUtilities;

/**
 * Construct a temporary Option object of a particular
 * type to determine what the appropriate invocation
 * of casm-utilities command should be.
 */

template <typename T>
std::string option_desc()
{
    T dumbopt;
    return dumbopt.tag();
}

void basic_help()
{
    std::cout<<"The available options for casm-utilities are:"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"    "<<casmUtilities::splash_launcher_name()<<std::endl;
    std::cout<<"    "<<casmUtilities::structure_launcher_name()<<std::endl;
    std::cout<<"    "<<casmUtilities::tester_launcher_name()<<std::endl;
    std::cout<<"    "<<casmUtilities::bedazzle_launcher_name()<<std::endl;
    std::cout<<"    "<<casmUtilities::boxify_launcher_name()<<std::endl;
    std::cout<<std::endl;

    return;
}

int main(int argc, char *argv[])
{
    if(argc==1)
    {
        basic_help();
        return 0;
    }

    if(argv[1]==casmUtilities::splash_launcher_name())
    {
        splash_utility_launch(argc, argv);
    }

    else if(argv[1]==casmUtilities::structure_launcher_name())
    {
        structure_utility_launch(argc, argv);
    }

    else if(argv[1]==casmUtilities::tester_launcher_name())
    {
        tester_utility_launch(argc, argv);
    }

    else if(argv[1]==casmUtilities::bedazzle_launcher_name())
    {
        bedazzle_utility_launch(argc, argv);
    }

    else if(argv[1]==casmUtilities::boxify_launcher_name())
    {
        boxify_utility_launch(argc, argv);
    }

    else
    {
        std::cout<<"Unrecognized option in casm-utilities!"<<std::endl;
        basic_help();
    }

    return 0;
}
