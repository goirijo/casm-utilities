#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/xtal/frankenstein.hpp"
#include "casmutils/xtal/structure_tools.hpp"

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace utilities
{

void frankenstack_initializer(po::options_description& frankenstack_desc)
{
    utilities::add_help_suboption(frankenstack_desc);
    utilities::add_desc_suboption(frankenstack_desc);
    utilities::add_output_suboption(frankenstack_desc);
    frankenstack_desc.add_options()("substructures,s", po::value<std::vector<fs::path>>()->multitoken()->required(),
                                    "POS.vasp like files you want to "
                                    "stack on top of each other.");
    frankenstack_desc.add_options()("number,n", po::value<int>(),
                                    "number of times to repeat the "
                                    "stacking if only one unit is given");
    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    Handler frankenstack_launch(argc, argv, frankenstack_initializer);

    if (frankenstack_launch.count("help"))
    {
        std::cout << frankenstack_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        frankenstack_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto sub_paths = frankenstack_launch.fetch<std::vector<fs::path>>("substructures");

    if (!sub_paths.size())
    {
        std::cerr << "Missing substructure path " << std::endl;
    }
    auto big_struc = rewrap::Structure::from_poscar(sub_paths[0]);
    if (sub_paths.size() == 1)
    {
        auto unit = rewrap::Structure::from_poscar(sub_paths[0]);
        int num_stacks = 1;
        if (frankenstack_launch.vm().count("number"))
        {
            num_stacks = frankenstack_launch.fetch<int>("number");
        }
        std::vector<rewrap::Structure> struc_vec;
        struc_vec.insert(struc_vec.end(), num_stacks, unit);
        big_struc = frankenstein::stack(struc_vec);
    }
    else
    {
        std::vector<rewrap::Structure> units;
        for (auto& item : sub_paths)
        {
            units.push_back(rewrap::Structure::from_poscar(item));
        }
        big_struc = frankenstein::stack(units);
    }
    std::cout << "                         .eeeeeeeee\n"
                 "	    		.$$$$$$$$P\n"
                 "		       .$$$$$$$$P\n"
                 "		      z$$$$$$$$P\n"
                 "		     z$$$$$$$$\n"
                 "		    z$$$$$$$$\n"
                 "		   d$$$$$$$$\n"
                 "		  d$$$$$$$$\n"
                 "		.d$$$$$$$P\n"
                 "	       .$$$$$$$$P\n"
                 "	      .$$$$$$$$$.........\n"
                 "	     .$$$$$$$$$$$$$$$$$$\n"
                 "	    z$$$$$$$$$$$$$$$$$P\n"
                 "	   -**********$$$$$$$P\n"
                 "		     d$$$$$$\n"
                 "		   .d$$$$$$\n"
                 "		  .$$$$$$P\n"
                 "		 z$$$$$$P\n"
                 "		d$$$$$$\n"
                 "	      .d$$$$$$\n"
                 "	     .$$$$$$$\n"
                 "	    z$$$$$$$beeeeee\n"
                 "	   d$$$$$$$$$$$$$*\n"
                 "	  ^^^^^^^^^$$$$$\n"
                 "		  d$$$*\n"
                 "		 d$$$\n"
                 "		d$$*\n"
                 "	       d$P\n"
                 "	     .$$\n"
                 "	    .$P\n"
                 "	   .$\n"
                 "	  .P\n"
                 "	 .\n"
                 "	/\n"
              << std::endl;
    std::cout << "IT'S ALIVE " << std::endl;
    if (frankenstack_launch.vm().count("output"))
    {
        casmutils::xtal::write_poscar(big_struc, frankenstack_launch.fetch<fs::path>("output"));
        return 0;
    }
    casmutils::xtal::print_poscar(big_struc, std::cout);
    return 0;
}
