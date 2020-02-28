#include <casmutils/definitions.hpp>
#include <casmutils/handlers.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace utilities
{

void primify_initializer(po::options_description& primify_desc)
{
    utilities::add_help_suboption(primify_desc);
    utilities::add_output_suboption(primify_desc);

    primify_desc.add_options()("superstructure,s", po::value<fs::path>()->required(),
                               "POS.vasp like file that you want to get the primitive structure for.");

    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    Handler primify_launch(argc, argv, primify_initializer);

    /* if(primify_launch.count("help") || primify_launch.argc()<2) */
    if (primify_launch.count("help"))
    {
        std::cout << primify_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        primify_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto super_path = primify_launch.fetch<fs::path>("superstructure");

    // Should all CASM calls be wrapped up?
    auto super_struc = rewrap::Structure::from_poscar(super_path);
    auto prim_struc = casmutils::xtal::make_primitive(super_struc);

    if (primify_launch.vm().count("output"))
    {
        auto out_path = primify_launch.fetch<fs::path>("output");
        casmutils::xtal::write_poscar(prim_struc, out_path);
    }

    else
    {
        casmutils::xtal::print_poscar(prim_struc, std::cout);
    }

    return 0;
}
