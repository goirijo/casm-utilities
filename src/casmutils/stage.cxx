#include "casmutils/stage.hpp"
#include <casm/crystallography/Structure.hh>
#include <casm/clex/ScelEnum.hh>
#include <casm/crystallography/SupercellEnumerator.hh>
#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/lattice.hpp"
#include <casm/crystallography/Niggli.hh>

namespace Simplicity
{

namespace
{
// surface area of a given lattice
double lattice_surface_area(const Rewrap::Lattice& lat)
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
double boxy_score(const Rewrap::Lattice& lat)
{
    // Less surface area per volume means more boxy
    // i.e. more volume per surface area means more boxy
    return std::abs(lat.vol()) / lattice_surface_area(lat);
}
} // namespace

// Finds the superstructure with the highest volume/surface_area
// Assuming that the input has structures of same volume
std::vector<Rewrap::Structure>::size_type boxiest_structure_index(const std::vector<Rewrap::Structure>& candidate_structures)
{
    //TODO: throw exception on empty vector
    double running_score = 0;
    std::vector<Rewrap::Structure>::size_type ix=0;
    std::vector<Rewrap::Structure>::size_type best_ix=ix;
    for (const auto& scel : candidate_structures)
    {
        double candidate_score = boxy_score(scel.lattice());
        if (candidate_score > running_score)
        {
            running_score = candidate_score;
            best_ix=ix;
        }
        ++ix;
    }
    return best_ix;
}

// Find the boxiest superstructure per volume for range of volumes
Rewrap::Structure make_boxiest_superstructure_of_volume(Rewrap::Structure& structure, int volume)
{
        std::vector<Rewrap::Structure> same_vol_scels = make_superstructures_of_volume(structure, volume);
        return same_vol_scels[boxiest_structure_index(same_vol_scels)];
}

std::vector<Rewrap::Structure> make_superstructures_of_volume(const Rewrap::Structure& structure, int volume)
{
    std::vector<Rewrap::Structure> all_superstructures;
    CASM::ScelEnumProps enum_props(volume, volume+1);
    CASM::SupercellEnumerator<CASM::Lattice> lat_enumerator(structure.lattice(), enum_props, CASM::TOL);

    for (const auto& lat : lat_enumerator)
    {
        Rewrap::Structure super = structure.create_superstruc(lat);
        Simplicity::make_niggli(&super);
        all_superstructures.emplace_back(std::move(super));
    }

    return all_superstructures;
}

} // namespace SuperBoxy
