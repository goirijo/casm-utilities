#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void superstructure_initializer(po::options_description& superstructure_desc)
{
    UtilityProgramOptions::add_help_suboption(superstructure_desc);
    UtilityProgramOptions::add_output_suboption(superstructure_desc);

    superstructure_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                               "POS.vasp like file that you want to get the super structure for.");
    superstructure_desc.add_options()("transf-matrix,t", po::value<fs::path>()->required(),
                                      "path to a file with transformation matrix.");

    return;
}
}

using namespace Utilities;

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
    Rewrap::fs::ifstream mat_file(transf_file_path);
    mat_file >> transf_mat;

    auto struc = Rewrap::Structure(structure_path);
    auto super_struc = Simplicity::make_super_structure(struc, transf_mat);

    //checks the output type and writes the super structure to a output stream
    if (superstructure_launch.vm().count("output"))
    {
        auto out_path = superstructure_launch.fetch<fs::path>("output");
        Simplicity::write_poscar(super_struc, out_path);
    }

    else
    {
        Simplicity::print_poscar(super_struc, std::cout);
    }

    return 0;
}
