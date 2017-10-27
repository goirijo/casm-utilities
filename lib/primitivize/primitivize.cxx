#include <iostream>
#include <fstream>

#include <casm/completer/Handlers.hh>
#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/Site.hh>
#include "lib/primitivize/primitivize.hpp"
#include "lib/completer/handlers.hpp"
#include "lib/launch/launchers.hpp"
#include "lib/alone/simplicity.hpp"


namespace casmUtilities
{
    std::string primitivize_launcher_name()
    {
        return "primitivize";
    }


    LaunchRuleList primitivize_initializer(po::options_description &primitivize_desc)
    {
        utilityProgramOptions::add_help_suboption(primitivize_desc);
        utilityProgramOptions::add_output_suboption(primitivize_desc);
        utilityProgramOptions::add_structure_suboption(primitivize_desc);

        /* primitivize_desc.add_options() */
        /*     ("check,c", */

        LaunchRuleList primitivize_rules;

        return primitivize_rules;
    }

    void primitivize_utility_launch(int argc, char *argv[])
    {
        Launcher primitivize_launch(argc, argv, primitivize_launcher_name(), primitivize_initializer);

        if (primitivize_launch.count("help") || primitivize_launch.argc()==2)
        {
            std::cout << primitivize_launch.utility().desc() << std::endl;
            return;
        }

        try
        {
            primitivize_launch.notify();
        }

        catch (boost::program_options::required_option &e)
        {
            std::cerr << e.what() << std::endl;
            return;
        }

        /* //Determine where to print output (if necessary) */
        /* std::ostream *out_stream=&std::cout; */
        /* std::ofstream outfile_stream; */ 
        /* if(primitivize_launch.count("output")) */
        /* { */
        /*     fs::path out_path=primitivize_launch.fetch<fs::path>("output"); */
        /*     outfile_stream.open(out_path.string()); */
        /*     out_stream=&outfile_stream; */
        /* } */
        auto target=primitivize_launch.fetch<CASM::fs::path>("output");
        std::ofstream target_stream(target.string());
        CASM::fs::path struc_path=primitivize_launch.fetch<fs::path>("structure");
        auto input_struc=simple::basic_structure_from_path(struc_path);

        CASM::BasicStructure<CASM::Site> prim_struc;
        if(!input_struc.is_primitive(prim_struc))
        {
            simple::structure_print(target_stream, prim_struc);
            std::cout<<"Primitive structre printed to "<<struc_path.string()<<std::endl;
        }

        else
        {
            std::cout<<"The given structure is already primitive!"<<std::endl;
        }



        return;
    }
}
