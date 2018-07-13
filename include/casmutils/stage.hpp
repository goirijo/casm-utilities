#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

#include <iostream>
#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>

/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction)
/// and splits it into n+1 structures where n is the size of slice_loc and the
/// entries in slice_loc
/// dictate where the slicing happens in increasing order. Entries in slice_loc
/// are fractional lengths of the c vector
std::vector<Rewrap::Structure> structure_slicer(const Rewrap::Structure &big_struc,
					std::vector<double> &slice_loc);

/// This function splits a structure in equally sized slices along the c -axis
/// (layering)
/// which is perpendicular to the ab plane.
std::vector<Rewrap::Structure> structure_slicer(const Rewrap::Structure &big_struc,
					int n_pieces);

/// This function takes a vector of structures with the same ab lattice vectors
/// and stacks them along the c direction. A corresponding vector of ab in plane
/// shift values relative to the first structure is given
Rewrap::Structure structure_stacker(
    std::vector<Rewrap::Structure> &sub_strucs,
    std::vector<std::pair<double, double>> &shift_vals);

/// This version of structure sets the same horizontal shift between each
/// consecutive layer
/// Layers are possibly different
Rewrap::Structure structure_stacker(std::vector<Rewrap::Structure> &sub_strucs,
			    const std::pair<double, double> &shift_value);

/// This version of structure sets the same layer as the stacking unit
/// Shifts are possibly different
Rewrap::Structure structure_stacker(
    Rewrap::Structure &unit, std::vector<std::pair<double, double>> &shift_vals);

/// This version of structure sets the same layer as the stacking unit and the
/// same shift between each layer n times
Rewrap::Structure structure_stacker(Rewrap::Structure &unit,
			    const std::pair<double, double> &shift_value,
			    int n_times);




#endif
