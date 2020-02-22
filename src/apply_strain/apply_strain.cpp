#include "casmutils/definitions.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/structure_tools.hpp"
#include <boost/filesystem/fstream.hpp> //pls no
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace Utilities
{

void applystrain_initializer(po::options_description& applystrain_desc)
{
    UtilityProgramOptions::add_help_suboption(applystrain_desc);
    UtilityProgramOptions::add_output_suboption(applystrain_desc);

    applystrain_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                   "POS.vasp like file that you want to apply strain to.");
    applystrain_desc.add_options()("mode,m", po::value<std::string>()->default_value("GL"),
                                   "Accepts strain convention as mode ('GL' [Green-Lagrange, default], 'EA' "
                                   "[Euler-Almansi], 'B' [Biot], or 'H' [Hencky])."
                                   " Also accepts 'F' [Deformation] as an argument to apply a deformation tensor");
    applystrain_desc.add_options()("tensor,t", po::value<fs::path>()->required(),
                                   "Path to a file with strain tensor."
                                   " Unrolled strain should be provided for GL, B, H, EA modes."
                                   " Ordered as E(0,0) E(1,1) E(2,2) E(1,2) E(0,2) E(0,1)."
                                   " Takes a 3X3 matrix for F (Deformation) mode.");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler applystrain_launch(argc, argv, applystrain_initializer);

    if (applystrain_launch.count("help"))
    {
        std::cout << applystrain_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        applystrain_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = applystrain_launch.fetch<fs::path>("structure");
    const auto strain_path = applystrain_launch.fetch<fs::path>("tensor");
    const auto mode = applystrain_launch.fetch<std::string>("mode");
    Rewrap::Structure strained_struc(struc_path);

    // check if the mode is a strain convention type or deformation mode
    // reads the input as a vector if its an unrolled strain in case of GL, B, H, EA modes else if its deformation mode
    // reads as a matrix
    if (mode == "F")
    {
        Eigen::Matrix3d deformation_tensor;
        boost::filesystem::ifstream mat_file(strain_path);
        mat_file >> deformation_tensor;
        Simplicity::apply_deformation(&strained_struc, deformation_tensor);
    }
    else
    {
        try
        {
            Eigen::VectorXd unrolled_strain(6);
            boost::filesystem::ifstream mat_file(strain_path);
            mat_file >> unrolled_strain;
            Simplicity::apply_strain(&strained_struc, unrolled_strain, mode);
        }

        catch (UtilExcept::UserInputMangle)
        {
            std::cerr << "Unrecognized --mode " << mode << std::endl;
            std::cerr << "options are GL(GREEN_LAGRANGE), B(BIOT), H(HENCKY), EA(EULER_ALMANSI), and F(Deformation)"
                      << std::endl;
            return 3;
        }
    }

    // checks the output type and writes the strained structure to a output stream
    if (applystrain_launch.vm().count("output"))
    {
        auto out_path = applystrain_launch.fetch<fs::path>("output");
        Simplicity::write_poscar(strained_struc, out_path);
    }

    else
    {
        Simplicity::print_poscar(strained_struc, std::cout);
    }

    return 0;
}
