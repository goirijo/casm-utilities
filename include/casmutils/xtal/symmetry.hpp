#ifndef XTAL_SYMMETRY_HH
#define XTAL_SYMMETRY_HH

#include <casmutils/sym/cartesian.hpp>
#include <vector>

namespace casmutils
{

namespace xtal
{
class Lattice;
class Structure;

/// Create the group of all symmetry operations that map the
/// given lattice onto itself.
std::vector<sym::CartOp> make_point_group(const Lattice& lat, double tol);

/// Create the group of symmetry operations that map the given
/// crystal structure onto itself.
std::vector<sym::CartOp> make_factor_group(const Structure& struc, double tol);

/// Modify the given Lattice such that it perfectly obeys the provided
/// symmetry group. Useful for reducing noise in lattice vectors.
Lattice symmetrize(const Lattice& noisy_lattice, const std::vector<sym::CartOp>& enforced_point_group);

// TODO: Does this even make sense?
/// Modify the given Structure such that it perfectly obeys the provided
/// symmetry group. Useful for reducing noise in lattice vectors and basis.
/// First, the group will by applied to the lattice, to even out the vectors.
/// The same operations will then be applied to the basis and the resulting values
/// will be averaged out.
/* Lattice symmetrize(const Lattice& noisy_lattice, const std::vector<sym::CartOp>& enforced_point_group); */

} // namespace xtal
} // namespace casmutils

#endif
