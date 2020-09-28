#ifndef FRANKENSTEIN_HH
#define FRANKENSTEIN_HH

#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <set>

namespace casmutils
{
/// The evil doctor has come to town the namespace frankenstein is
/// for doing all the things the doctor did to that poor monster.
/// Snipping structures into partial components and Fastening them back together
/// to create whatever atrocity you seek
namespace frankenstein
{

/// This function takes a structure (where the slice direction is perpendicular to ab plane)
/// and splits it into two structures.
/// slice_loc is in fractional length of the c vector
/// the first item in the pair is below the slice location to the 0 boundary,
/// the second is above the slice location to the 1 boundary
std::pair<xtal::Structure, xtal::Structure> slice(const xtal::Structure& big_struc, double slice_loc, double tol);

/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction)
/// and splits it into n+1 structures where n is the size of slice_locs and the
/// entries in slice_locs
/// dictate where the slicing happens in increasing order. Entries in slice_loc
/// are fractional lengths of the c vector
std::vector<xtal::Structure> multi_slice(const xtal::Structure& big_struc,
                                         const std::set<double>& slice_locs, // consider using std::set
                                         double tol);

/// This function splits a structure in equally sized slices along the c -axis
/// (layering)
/// which is perpendicular to the ab plane.
std::vector<xtal::Structure> uniformly_slice(const xtal::Structure& big_struc, int n_pieces);

/// This function takes a structure and reduces the lattice boundaries to the
/// closest atoms to each of the specified boundaries dictated by the vector
/// dirs. dirs is a vector of bools that indicate whether or not to shrink
/// along the a, b, and c direction respectively.
/// An amount of extra space is left, specified by the padding value, so that
/// atoms at the boundaries don't overlap in the periodic images.
xtal::Structure vacuum_pack(const xtal::Structure& big_struc, std::array<bool, 3>& dirs, double padding);

/// This function takes a structure and increases the lattice boundaries by
/// shift value. Shift value is the extension in cartesian length along
/// along the a, b, and c direction respectively.
xtal::Structure inflate(const xtal::Structure& struc, const std::array<double, 3>& padding);

/// This function takes a structures and shifts the origin by shift val
/// shift val is in fractional coordinates of the lattice
/* void shift_coords_by(xtal::Structure* struc, const Eigen::Vector3d& shift_val); */

/// Stack a series of structres by adjusting the ab-vectors to match, and then concatenatig
/// along the c axis
xtal::Structure stack(const std::vector<xtal::Structure>& sub_strucs);

/// Translate the given basis by the specified cartesian value
std::vector<xtal::Site> translate_basis(const std::vector<xtal::Site>& basis, const Eigen::Vector3d& shift);

/// Translate the entire basis of the structure by the specified amount
xtal::Structure translate_basis(const xtal::Structure& struc, const Eigen::Vector3d& shift);
} // namespace frankenstein
} // namespace casmutils
#endif
