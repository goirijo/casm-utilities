#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>
#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"

/// The evil doctor has come to town the namespace Frankenstein is
/// for doing all the things the doctor did to that poor monster.
/// Snipping structures into partial components and Fastening them back together
/// to create whatever atrocity you seek
namespace Frankenstein {

/// This function takes a structure (where lattice vector c is perpendicular to ab plane and
/// is layering direction) and splits it into two structures.
/// slice_loc is in fractional length of the c vector
/// the first item in the pair is below the slice location to the 0 boundary,
/// the second is above the slice location to the 1 boundary
std::pair<Rewrap::Structure, Rewrap::Structure> slice(
    const Rewrap::Structure &big_struc, double slice_loc, double tol);

/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction)
/// and splits it into n+1 structures where n is the size of slice_locs and the
/// entries in slice_locs
/// dictate where the slicing happens in increasing order. Entries in slice_loc
/// are fractional lengths of the c vector
std::vector<Rewrap::Structure> multi_slice(const Rewrap::Structure &big_struc,
					   const std::set<double> &slice_locs,    //consider using std::set
					   double tol);

/// This function splits a structure in equally sized slices along the c -axis
/// (layering)
/// which is perpendicular to the ab plane.
std::vector<Rewrap::Structure> uniformly_slice(const Rewrap::Structure &big_struc,
					   int n_pieces);

/// This function takes a vector of structures with the same ab lattice vectors
/// and stacks them along the c direction.
Rewrap::Structure stack(
    const std::vector<Rewrap::Structure> &sub_strucs);

/// This function takes a structure and reduces the lattice boundaries to the
/// closest atoms to each of the specified boundaries dictated by the vector
/// dirs. dirs is a vector of bools that indicate whether or not to shrink
/// along the a, b, and c direction respectively.
/// An amount of extra space is left, specified by the padding value, so that
/// atoms at the boundaries don't overlap in the periodic images.
Rewrap::Structure vacuum_pack(const Rewrap::Structure &big_struc,
			      std::array<bool,3> &dirs, double padding);

/// This function takes a structure and increases the lattice boundaries by
/// shift value. Shift value is the extension in cartesian length along
/// along the a, b, and c direction respectively.
Rewrap::Structure inflate(const Rewrap::Structure &struc,
			  const std::array<double,3> &padding);

/// This function takes a structures and shifts the origin by shift val
/// shift val is in fractional coordinates of the lattice
void shift_coords_by(Rewrap::Structure *struc,
				const Eigen::Vector3d &shift_val);
}

#endif
