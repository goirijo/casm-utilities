#ifndef MISC_HH
#define MISC_HH

#include <casm/external/Eigen/Core>
#include <casm/misc/CASM_Eigen_math.hh>
#include <casm/misc/type_traits.hh>
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
using CASM::almost_equal;
using CASM::almost_zero;

template <typename ComparatorType_f, typename CompareType, typename... Args>
bool is_equal(const CompareType& reference, const CompareType& other, const Args&... functor_params)
{
    ComparatorType_f reference_equals(functor_params...);
    return reference_equals(reference, other);
}

/// Use this class to convert a BinaryComparator_f functor to UnaryComparator_f functor useful in std::find_if
/// For example if you have a "LatticeEquals_f" BinaryComparator and you want a UnaryComparator_f out of it:
/// UnaryComparator_f<LatticeEquals_f> lattice_unary_compare(ref_lattice, tol);
template <typename BinaryComparator_f> class UnaryComparator_f
{
public:
    using CompareType = notstd::first_argument_type<BinaryComparator_f>;

    template <typename... Args>
    UnaryComparator_f(const CompareType& ref_value, const Args&... functor_params)
        : m_ref_value(ref_value), m_compare_method(functor_params...)
    {
    }

    bool operator()(const CompareType& other_value) { return m_compare_method(m_ref_value, other_value); }

private:
    CompareType m_ref_value;
    BinaryComparator_f m_compare_method;
};

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
} // namespace extend

namespace io
{
Eigen::IOFormat coord_format();
}

#endif
