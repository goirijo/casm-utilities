#include "casmutils/stage.hpp"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include <boost/filesystem.hpp>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <set>

namespace SuperBoxy
{

std::vector<Rewrap::Lattice> make_superlattices_of_size(Rewrap::Structure& tile, int volume);

std::vector<Rewrap::Structure> make_supercells(Rewrap::Structure& structure, int min_vol, int max_vol)
{
    std::vector<Rewrap::Structure> all_supercells;
    int max_vol_incl = max_vol + 1;
    CASM::ScelEnumProps enum_props(min_vol, max_vol_incl);
    CASM::SupercellEnumerator<CASM::Lattice> lat_enumerator(structure.lattice(), enum_props);

    for(const auto& lat : lat_enumerator)
    {
        //TODO: Make niggli before creating structure
        auto super = Simplicity::make_niggli(structure.create_superstruc(lat));
        all_supercells.push_back(super);
    }

    return all_supercells;
}

// surface area from lattice borrowed from John Goiri
double lattice_surface_area(const CASM::Lattice& lat)
{
    Eigen::Vector3d a = lat[0];
    Eigen::Vector3d b = lat[1];
    Eigen::Vector3d c = lat[2];

    double ab = a.cross(b).norm();
    double bc = b.cross(c).norm();
    double ca = c.cross(b).norm();

    return std::abs(ab) + std::abs(bc) + std::abs(ca);
}

// score for determining level of boxiness, borrowed from John Goiri
double boxy_score(const CASM::Lattice& lat)
{
    // Less surface area per volume means more boxy
    // i.e. more volume per surface area means more boxy
    return std::abs(lat.vol()) / lattice_surface_area(lat);
}

// Finds the supercell with the highest volume/surface_area
// Assuming that the input has structures of same volume
Rewrap::Structure most_boxy(std::vector<Rewrap::Structure>& supercells)
{
    double running_score = 0;
    Rewrap::Structure boxiest_scel = supercells[0];
    for (const auto& scel : supercells)
    {
        double candidate_score = boxy_score(scel.lattice());
        if (candidate_score > running_score)
        {
            running_score = candidate_score;
            boxiest_scel = scel;
        }
    }

    return boxiest_scel;
}

// Find the boxiest supercell per volume for range of volumes
std::vector<Rewrap::Structure> make_boxy_supercells(Rewrap::Structure& structure, int min_vol, int max_vol)
{
    std::vector<Rewrap::Structure> boxy_supercells;
    int max_vol_incl = max_vol + 1;
    for (int a = min_vol; a < max_vol_incl; ++a)
    {
        std::vector<Rewrap::Structure> same_vol_scels = make_supercells(structure, a, a);
        Rewrap::Structure same_vol_boxy = most_boxy(same_vol_scels);
        boxy_supercells.push_back(same_vol_boxy);
    }
    return boxy_supercells;
}
} // namespace SuperBoxy
