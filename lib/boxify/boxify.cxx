#include <iostream>
#include <string>

#include "lib/completer/handlers.hpp"
#include "lib/definitions.hpp"
#include "lib/launch/launchers.hpp"
#include "lib/alone/simplicity.hpp"

#include <casm/completer/Complete.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/crystallography/SupercellEnumerator.hh>

namespace casmUtilities
{
namespace boxifyImpl
{
double lattice_surface_area(const CASM::Lattice &lat)
{
    Eigen::Vector3d a = lat[0];
    Eigen::Vector3d b = lat[1];
    Eigen::Vector3d c = lat[2];

    double ab = a.cross(b).norm();
    double bc = b.cross(c).norm();
    double ca = c.cross(b).norm();

    return std::abs(ab) + std::abs(bc) + std::abs(ca);
}

double boxy_score(const CASM::Lattice &lat)
{
    // Less surface area per volume means more boxy
    // i.e. more volume per surface area means more boxy
    return std::abs(lat.vol()) / lattice_surface_area(lat);
}

CASM::Index boxiest_supercell_index(const CASM::PrimClex &pclex)
{
    double running_score=0;
    CASM::Index boxiest_ix;
    for(CASM::Index i=1; i<pclex.get_supercell_list().size(); ++i)
    {
        const auto& prim_lat=pclex.get_prim().lattice();
        const auto& transf_mat=pclex.get_supercell(i).get_transf_mat();
        CASM::Lattice candidate_lat=CASM::make_supercell(prim_lat, transf_mat);

        double candidate_score=boxy_score(candidate_lat);
        if(candidate_score>running_score)
        {
            running_score=candidate_score;
            boxiest_ix=i;
        }
    }

    return boxiest_ix;
}

}

std::string boxify_launcher_name() { return "boxify"; }

LaunchRuleList boxify_initializer(po::options_description &boxify_desc)
{
    utilityProgramOptions::add_help_suboption(boxify_desc);

    boxify_desc.add_options()("structure,s",
                              po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path())->required(),
                              "POS.vasp like file that you want to create a boxy supercell for.")(
        "volume,v", po::value<int>()->required(), "Specifies the size of the resulting supercell.")(
        "output,o", po::value<std::string>()->required(), "Output file for structure");

    LaunchRuleList boxify_rules;

    return boxify_rules;
}

void boxify_utility_launch(int argc, char *argv[])
{
    Launcher boxify_launch(argc, argv, boxify_launcher_name(), boxify_initializer);

    if (boxify_launch.count("help"))
    {
        std::cout << boxify_launch.utility().desc() << std::endl;
        return;
    }

    try
    {
        boxify_launch.notify();
    }

    catch (boost::program_options::required_option &e)
    {
        std::cerr << e.what() << std::endl;
        return;
    }

    //Read structure and initialize primclex with it as prim
    CASM::fs::path struc_path=boxify_launch.fetch<fs::path>("structure");
    auto base_struc=simple::basic_structure_from_path(struc_path);
    CASM::Structure base_struc_recast(base_struc);
    CASM::PrimClex pclex(base_struc_recast);

    //Determine the desired volume and create supercell enumerator
    int volume=boxify_launch.fetch<int>("volume");
    CASM::ScelEnumProps enumeration(volume,volume+1);
    pclex.generate_supercells(enumeration);

    //Determine the boxies supercell of them all and get the corresponding superstructure
    auto boxy_struc=pclex.get_supercell(boxifyImpl::boxiest_supercell_index(pclex)).superstructure();

    //Print the result out!
    std::string target=boxify_launch.fetch<std::string>("output");
    std::ofstream target_stream(target);
    /* simple::structure_print(std::cout, boxy_struc); */
    simple::structure_print(target_stream, boxy_struc);

    return;
}
}
