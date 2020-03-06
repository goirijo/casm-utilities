#ifndef MISC_HH
#define MISC_HH

#include <casm/external/Eigen/Core>
#include <string>
#include <vector>

namespace CASM
{
namespace xtal
{
class Site;
class Coordinate;
class Structure;
} // namespace xtal
class PrimClex;
} // namespace CASM

namespace casmutils
{

template <typename ComparatorType_f, typename CompareType, typename... Args>
bool is_equal(const CompareType& reference, const CompareType& other, const Args&... functor_params)
{
    ComparatorType_f reference_equals(reference, functor_params...);
    return reference_equals(other);
}
} // namespace casmutils

/**
 * This namespace is reserved for convenience functions
 * that reduce boilerplate code within library functions
 * of casm-utilities (e.g. simplicity).
 * Usage of anything within this namespace should not leak
 * outside of any implementation, these are not utility
 * library functions, they're just here for convenience.
 */
namespace extend
{
/// Make a PrimClex but shut up about it
CASM::PrimClex quiet_primclex(CASM::xtal::Structure& prim);

/// Alternate construction of a site with only occupant degrees of freedom, and atom as occuapnt (no molecule)
CASM::xtal::Site atomic_site(const CASM::xtal::Coordinate& coord, const std::vector<std::string>& allowed_species);
} // namespace extend

namespace io
{
Eigen::IOFormat coord_format();
}

#endif
