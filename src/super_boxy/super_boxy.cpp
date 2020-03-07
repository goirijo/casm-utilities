#include <casmutils/definitions.hpp>
#include <casmutils/handlers.hpp>
#include <casmutils/stage.hpp>
#include <casmutils/xtal/structure_tools.hpp>

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace utilities
{

void super_boxy_initializer(po::options_description& super_boxy_desc)
{
    utilities::add_help_suboption(super_boxy_desc);
    utilities::add_output_suboption(super_boxy_desc);
    super_boxy_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                  "POS.vasp like file that you want to get the boxy supercell for.");
    super_boxy_desc.add_options()("volume,v", po::value<int>()->required(),
                                  "Volume of the boxy superstructure, relative to the input structure");

    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    using namespace casmutils;
    Handler super_boxy_launch(argc, argv, super_boxy_initializer);

    if (super_boxy_launch.count("help"))
    {
        std::cout << super_boxy_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        super_boxy_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    // initialize structure from POSCAR file
    // xtal::Structure in_struc("/home/julija/programming/casm-utilities/test/POSCAR_1");
    auto in_struc_path = super_boxy_launch.fetch<fs::path>("structure");
    auto in_struc = xtal::Structure::from_poscar(in_struc_path);
    auto in_vol = super_boxy_launch.fetch<int>("volume");

    auto boxy_struc = casmutils::xtal::make_boxiest_superstructure_of_volume(in_struc, in_vol);

    if (super_boxy_launch.vm().count("output"))
    {
        auto out_path = super_boxy_launch.fetch<fs::path>("output");
        casmutils::xtal::write_poscar(boxy_struc, out_path);
    }

    else
    {
        casmutils::xtal::print_poscar(boxy_struc, std::cout);
    }

    return 0;
}
