#include "casmutils/definitions.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void super_boxy_initializer(po::options_description& super_boxy_desc)
{
    UtilityProgramOptions::add_help_suboption(super_boxy_desc);
    UtilityProgramOptions::add_output_suboption(super_boxy_desc);

    //super_boxy_desc.add_options()("superstructure,s", po::value<fs::path>()->required(),
//                               "POS.vasp like file that you want to get the boxy supercell for.");

    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
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
    //Rewrap::Structure my_structure("/home/julija/programming/casm-utilities/test/POSCAR_1");
    boost::filesystem::path my_path = "/home/julija/programming/casm-utilities/test/POSCAR_1";
    auto my_structure = Rewrap::Structure(my_path);

   //find supercells
    std::vector<Rewrap::Structure> all_supercells = SuperBoxy::make_supercells(my_structure, 1, 2); 
    for(const auto &super: all_supercells)
    { 
        Simplicity::print_poscar(super, std::cout);
    }

    std::cout << " ------------------------------------------------------- " << std::endl;
    std::cout << " Most boxy supercells " << std::endl;
    std::vector<Rewrap::Structure> boxy_supercells = SuperBoxy::make_boxy_supercells(my_structure, 1, 3); 
    for(const auto &super: boxy_supercells)
    { 
        Simplicity::print_poscar(super, std::cout);
    }


/*
    auto super_path = primify_launch.fetch<fs::path>("superstructure");

    // Should all CASM calls be wrapped up?
    auto super_struc = CASM::Structure(super_path);
    auto prim_struc = Simplicity::make_primitive(super_struc);

    if (primify_launch.vm().count("output"))
    {
        auto out_path = primify_launch.fetch<fs::path>("output");
        Simplicity::write_poscar(prim_struc, out_path);
    }

    else
    {
        Simplicity::print_poscar(prim_struc, std::cout);
    }
/**/
    return 0;
}
