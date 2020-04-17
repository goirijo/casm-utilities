#include <casmutils/definitions.hpp>
#include <casmutils/handlers.hpp>
#include <casmutils/xtal/frankenstein.hpp>
#include <casmutils/xtal/structure_tools.hpp>

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace utilities
{

void inflate_initializer(po::options_description& inflate_desc)
{
    utilities::add_help_suboption(inflate_desc);
    utilities::add_desc_suboption(inflate_desc);
    utilities::add_output_suboption(inflate_desc);
    inflate_desc.add_options()("structure,s",
                               po::value<fs::path>()->required(),
                               "POS.vasp like file you want to "
                               "increase the boundaries of");
    inflate_desc.add_options()("length,l",
                               po::value<std::vector<double>>()->multitoken()->required(),
                               "amount that each lattice vector will grow (three values "
                               "corresponding to a, b, and c)");
    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    using namespace casmutils;
    Handler inflate_launch(argc, argv, inflate_initializer);

    if (inflate_launch.count("help"))
    {
        std::cout << inflate_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        inflate_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto raw_lengths = inflate_launch.fetch<std::vector<double>>("length");

    if (raw_lengths.size() != 3)
    {
        std::cerr << "You must give exactly 3 --lengths for inflation, one value for each lattice direction."
                  << std::endl;
        return 3;
    }

    auto struc_path = inflate_launch.fetch<fs::path>("structure");
    auto struc = xtal::Structure::from_poscar(struc_path);
    auto out_struc = struc;
    std::array<double, 3> lengths{raw_lengths[0], raw_lengths[1], raw_lengths[2]};
    out_struc = frankenstein::inflate(struc, lengths);
    if (inflate_launch.vm().count("output"))
    {
        casmutils::xtal::write_poscar(out_struc, inflate_launch.fetch<fs::path>("output"));
        return 0;
    }
    casmutils::xtal::print_poscar(out_struc, std::cout);
    return 0;
}
