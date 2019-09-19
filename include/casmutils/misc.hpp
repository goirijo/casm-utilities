#ifndef MISC_HH
#define MISC_HH

#include <casm/external/Eigen/Core>
#include <string>
#include <vector>

namespace CASM
{
class Site;
class Coordinate;
class PrimClex;
class Structure;
} // namespace CASM

/**
 * This namespace is reserved for convenience functions
 * that reduce boilerplate code within library functions
 * of casm-utilities (e.g. Simplicity).
 * Usage of anything within this namespace should not leak
 * outside of any implementation, these are not utility
 * library functions, they're just here for convenience.
 */

namespace Extend
{
/// Make a PrimClex but shut up about it
CASM::PrimClex quiet_primclex(CASM::Structure& prim);

/// Alternate construction of a site with only occupant degrees of freedom, and atom as occuapnt (no molecule)
CASM::Site atomic_site(const CASM::Coordinate& coord, const std::vector<std::string>& allowed_species);
} // namespace Extend

namespace IO
{
Eigen::IOFormat coord_format();
}

#endif
