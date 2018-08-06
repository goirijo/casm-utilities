#include "casmutils/definitions.hpp"
#include "casmutils/frankenstein.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void coord_shift_initializer(po::options_description& coord_shift_desc)
{
    UtilityProgramOptions::add_help_suboption(coord_shift_desc);
    UtilityProgramOptions::add_desc_suboption(coord_shift_desc);
    UtilityProgramOptions::add_output_suboption(coord_shift_desc);
    coord_shift_desc.add_options()("structure", po::value<fs::path>()->required(),
                                   "POS.vasp like file you want to "
                                   "shift all the coordinates of");
    coord_shift_desc.add_options()("shift", po::value<std::vector<double>>()->multitoken()->required(),
                                   "shift value that will be added to all coordinates (negative "
                                   "origin shift) units are fractional");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler coord_shift_launch(argc, argv, coord_shift_initializer);

    if (coord_shift_launch.count("help"))
    {
        std::cout << coord_shift_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        coord_shift_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = coord_shift_launch.fetch<fs::path>("structure");
    auto struc = Rewrap::Structure(struc_path);
    auto out_struc = struc;
    auto vec = coord_shift_launch.fetch<std::vector<double>>("shift");
    Frankenstein::shift_coords_by(&out_struc, Eigen::Map<Eigen::Vector3d>(&vec[0]));
    if (coord_shift_launch.vm().count("output"))
    {
        Simplicity::write_poscar(out_struc, coord_shift_launch.fetch<fs::path>("output"));
        return 0;
    }
    Simplicity::print_poscar(out_struc, std::cout);
    return 0;
}
