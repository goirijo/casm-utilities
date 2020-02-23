#include "casmutils/definitions.hpp"
#include "casmutils/xtal/frankenstein.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/xtal/structure_tools.hpp"

#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace utilities
{

void frankenslice_initializer(po::options_description& frankenslice_desc)
{
    utilities::add_help_suboption(frankenslice_desc);
    utilities::add_desc_suboption(frankenslice_desc);
    utilities::add_output_suboption(frankenslice_desc);
    frankenslice_desc.add_options()("superstructure,s", po::value<fs::path>()->required(),
                                    "POS.vasp like file that you want to "
                                    "get the primitive structure for.");
    frankenslice_desc.add_options()("slice-locations,x", po::value<std::vector<double>>()->multitoken(),
                                    "a set of locations to cut the superstructure along, units are "
                                    "in fractional lengths of the c-axis");
    frankenslice_desc.add_options()("number,n", po::value<int>(), "number of equally sized pieces along c-axis");

    return;
}
} // namespace utilities

using namespace utilities;

int main(int argc, char* argv[])
{
    double tol = 1e-5;
    Handler frankenslice_launch(argc, argv, frankenslice_initializer);

    if (frankenslice_launch.count("help"))
    {
        std::cout << frankenslice_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        frankenslice_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto super_path = frankenslice_launch.fetch<fs::path>("superstructure");

    auto super_struc = rewrap::Structure::from_poscar(super_path);
    std::vector<rewrap::Structure> snippets;
    if (frankenslice_launch.vm().count("slice-locations"))
    {
        auto raw_slice_locs = frankenslice_launch.fetch<std::vector<double>>("slice-locations");
        std::set<double> slice_locs(raw_slice_locs.begin(), raw_slice_locs.end());
        snippets = frankenstein::multi_slice(super_struc, slice_locs, tol);
    }
    else if (frankenslice_launch.vm().count("number"))
    {
        auto num_slices = frankenslice_launch.fetch<int>("number");
        snippets = frankenstein::uniformly_slice(super_struc, num_slices);
    }
    else
    {
        std::cerr << "Neither vector or number option was given to "
                     "frankenslice"
                  << std::endl;
        return 3;
    }
    if (frankenslice_launch.vm().count("output"))
    {
        auto out_path = frankenslice_launch.fetch<fs::path>("output");
        int count = 0;
        for (auto& item : snippets)
        {
            // TODO: what if directory doesn't exist?
            std::ostringstream ostr;
            ostr << std::setfill('0') << std::setw(2) << count;
            simplicity::write_poscar(item, out_path / rewrap::fs::path("slice" + ostr.str() + "POSCAR"));
            count++;
        }
    }

    else
    {
        int count = 0;
        for (auto& item : snippets)
        {
            std::cout << "slice " << count << std::endl;
            simplicity::print_poscar(item, std::cout);
            count++;
        }
    }

    return 0;
}
