#include "casmutils/definitions.hpp"
#include "casmutils/frankenstein.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/structure_tools.hpp"

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace Utilities
{

void vacuumpack_initializer(po::options_description& vacuumpack_desc)
{
    UtilityProgramOptions::add_help_suboption(vacuumpack_desc);
    UtilityProgramOptions::add_desc_suboption(vacuumpack_desc);
    UtilityProgramOptions::add_output_suboption(vacuumpack_desc);
    vacuumpack_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                  "POS.vasp like file you want to "
                                  "shrink the boundaries of");
    vacuumpack_desc.add_options()("dirs,d", po::value<std::string>()->required(),
                                  "directions that shrinkage is allowed to happen in. any "
                                  "combination of a, b, and c");
    vacuumpack_desc.add_options()("padding,p", po::value<double>()->required(),
                                  "Add an amount of extra space around the border of the atom enclosure");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler vacuumpack_launch(argc, argv, vacuumpack_initializer);

    if (vacuumpack_launch.count("help"))
    {
        std::cout << vacuumpack_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        vacuumpack_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = vacuumpack_launch.fetch<fs::path>("structure");
    auto struc = Rewrap::Structure(struc_path);
    auto out_struc = struc;
    auto dirs = vacuumpack_launch.fetch<std::string>("dirs");
    auto padding = vacuumpack_launch.fetch<double>("padding");
    std::array<bool, 3> allowed_dirs;
    allowed_dirs[0] = (dirs.find("a") != std::string::npos);
    allowed_dirs[1] = (dirs.find("b") != std::string::npos);
    allowed_dirs[2] = (dirs.find("c") != std::string::npos);
    out_struc = Frankenstein::vacuum_pack(struc, allowed_dirs, padding);
    if (vacuumpack_launch.vm().count("output"))
    {
        Simplicity::write_poscar(out_struc, vacuumpack_launch.fetch<fs::path>("output"));
        return 0;
    }
    Simplicity::print_poscar(out_struc, std::cout);
    return 0;
}
