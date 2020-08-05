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
