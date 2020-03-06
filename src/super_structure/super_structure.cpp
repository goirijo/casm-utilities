#include "casmutils/exceptions.hpp"
#include <casmutils/definitions.hpp>
#include <casmutils/handlers.hpp>
#include <casmutils/xtal/structure_tools.hpp>

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace utilities
{

void superstructure_initializer(po::options_description& superstructure_desc)
{
    utilities::add_help_suboption(superstructure_desc);
    utilities::add_output_suboption(superstructure_desc);

    superstructure_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                      "POS.vasp like file that you want to get the super structure for.");
    superstructure_desc.add_options()("transf-matrix,t", po::value<fs::path>()->required(),
                                      "path to a file with transformation matrix.");

    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    Handler superstructure_launch(argc, argv, superstructure_initializer);

    if (superstructure_launch.count("help"))
    {
        std::cout << superstructure_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        superstructure_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto structure_path = superstructure_launch.fetch<fs::path>("structure");
    auto transf_file_path = superstructure_launch.fetch<fs::path>("transf-matrix");
    // read the matrix from the file into an eigen matrix
    Eigen::Matrix3i transf_mat;
    std::ifstream mat_file(
        transf_file_path); // If you came here looking for something weird going on, this is probably the cuplit
    mat_file >> transf_mat;

    throw except::NotImplemented();

    auto struc = rewrap::Structure::from_poscar(structure_path);
    auto super_struc = casmutils::xtal::make_super_structure(struc, transf_mat);

    // checks the output type and writes the super structure to a output stream
    if (superstructure_launch.vm().count("output"))
    {
        auto out_path = superstructure_launch.fetch<fs::path>("output");
        casmutils::xtal::write_poscar(super_struc, out_path);
    }

    else
    {
        casmutils::xtal::print_poscar(super_struc, std::cout);
    }

    return 0;
}
