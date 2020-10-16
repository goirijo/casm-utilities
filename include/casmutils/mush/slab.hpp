#ifndef MUSH_SLAB_HH
#define MUSH_SLAB_HH

#include "casmutils/xtal/site.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>
#include <vector>

namespace casmutils
{
namespace mush
{
/// Returns orhtogonal unit vectors oriented such that the point along the
/// a vector, the ab plane normal, and whatever is perpendicular to that (column vectors).
Eigen::Matrix3d slab_unit_vectors(const xtal::Lattice& slab);

/// Applies the given transformation the the *column* vector matrix representation
/// of the lattice
xtal::Lattice make_transformed_lattice(const xtal::Lattice& lat, const Eigen::Matrix3d& transform);

/// Returns the same lattice, but rotated such that the a vector points along the
/// Cartesian x direction, and the b vector is parallel to the xy plane.
xtal::Lattice make_aligned(const xtal::Lattice& struc);

/// Returns the same lattice, but rotated such that the a vector points along the
/// Cartesian x direction, and the b vector is parallel to the xy plane.
xtal::Structure make_aligned(xtal::Structure struc);
void make_aligned(xtal::Structure* struc);

/// Attempts to make the c vector more perpendicular to the ab-plane by applying a
///"within" operation to bring it's projection into the parallelogram spanned by the
/// a and b vectors.
xtal::Lattice orthogonalize_c_vector(const xtal::Lattice& lat);
xtal::Structure orthogonalize_c_vector(const xtal::Structure& lat);

/// Given a slab, where the surface plane has already been exposed to the a-b vectors, create a
/// superstructure by stacking units along the c direction
xtal::Structure make_stacked_slab(const xtal::Structure& slab_unit, int stacks);

/// Tranlsate the basis in the given structure such that the basis at the specified
/// index ends up at the origin
xtal::Structure make_floored_structure(const xtal::Structure& shiftable_struc, int floor_atom_ix);

// TODO:
/// Given a stacked slab lattice, reduce angles as much as possible in the c direction but still
/// keep the a-b vectors invariant
/* xtal::Lattice make_reduced_angles(const xtal::Lattice& primitive_unit, */
/*                                       const xtal::Lattice& slab_unit_lattice); */

} // namespace mush
} // namespace casmutils

#endif
