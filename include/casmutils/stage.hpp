#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <iostream>
#include "casmutils/structure.hpp"

namespace SuperBoxy
{
// Given a structure, find all the supercells between volumes min_vol and max_vol
        std::vector<Rewrap::Structure> make_supercells(Rewrap::Structure &structure, int min_vol, int max_vol);

// Find surface area given a lattice
    double lattice_surface_area(const CASM::Lattice &lat);

// volume to surface_area ratio
    double boxy_score(const CASM::Lattice &lat);

// Find the supercell with the highest volume/surface_area ratio of the ones given
    Rewrap::Structure most_boxy(std::vector<Rewrap::Structure> &supercells);

// Find the most boxy supercell at each volume
    std::vector<Rewrap::Structure> make_boxy_supercells(Rewrap::Structure &structure, int min_vol, int max_vol);
}
#endif
