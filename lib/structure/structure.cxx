#include <iostream>
#include <fstream>

#include <casm/completer/Handlers.hh>
#include "lib/structure/structure.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"


namespace casmUtilities
{
    std::string structure_launcher_name()
    {
        return "structure";
    }


    LaunchRuleList structure_initializer(po::options_description &structure_desc)
    {
        utilityProgramOptions::add_help_suboption(structure_desc);
        utilityProgramOptions::add_output_suboption(structure_desc);

        structure_desc.add_options()
        ("make-primitive,p",po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path()),"Convert the given POSCAR file into a primitive cell");

        LaunchRuleList structure_rules;
        structure_rules.add_any_inclusion("output",std::vector<std::string> {"make-primitive"});

        return structure_rules;
    }

    void structure_utility_launch(int argc, char *argv[])
    {
        Launcher structure_launch(argc, argv, structure_launcher_name(), casmUtilities::structure_initializer);

        //Determine where to print output (if necessary)
        std::ostream *out_stream=&std::cout;
        std::ofstream outfile_stream; 
        if(structure_launch.count("output"))
        {
            fs::path out_path=structure_launch.fetch<fs::path>("output");
            outfile_stream.open(out_path.string());
            out_stream=&outfile_stream;
        }

        if(structure_launch.count("make-primitive"))
        {
            *out_stream<<"Made primitive..."<<std::endl;
        }

        return;
    }
}
