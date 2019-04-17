#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"

namespace Simplicity
{
/// Given a structure, find all the superstructures between volumes min_vol and max_vol
std::vector<Rewrap::Structure> make_superstructures_of_volume(const Rewrap::Structure& structure, int volume);

/// Find the index of the superstructure with the highest volume/surface_area ratio of the ones given
std::vector<Rewrap::Structure>::size_type boxiest_structure_index(const std::vector<Rewrap::Structure>& candidate_structures);
/* const Rewrap::Structure& boxiest_structure(const std::vector<Rewrap::Structure>& candidate_structures); */

/// Find the most boxy superstructure at each volume
Rewrap::Structure make_boxiest_superstructure_of_volume(Rewrap::Structure& structure, int volume);
} // namespace SuperBoxy
#endif
