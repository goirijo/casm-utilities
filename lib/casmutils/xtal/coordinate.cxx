#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/lattice.hpp"

namespace casmutils
{
namespace xtal
{
Eigen::Vector3d fractional_to_cartesian(const Eigen::Vector3d& fractional_coordinates, const Lattice& lat)
{
    return CASM::xtal::Coordinate(fractional_coordinates, lat.__get(), CASM::FRAC).cart();
}

Eigen::Vector3d cartesian_to_fractional(const Eigen::Vector3d& cartesian_coord, const Lattice& lat)
{
    return CASM::xtal::Coordinate(cartesian_coord, lat.__get(), CASM::CART).frac();
}

Eigen::Vector3d bring_within_lattice(const Eigen::Vector3d& cartesian_coord, const Lattice& lat)
{
    CASM::xtal::Coordinate casm_coord(cartesian_coord, lat.__get(), CASM::CART);
    casm_coord.within();
    return casm_coord.cart();
}

Eigen::Vector3d bring_within_wigner_seitz(const Eigen::Vector3d& cartesian_coord, const Lattice& lat)
{
    CASM::xtal::Coordinate casm_coord(cartesian_coord, lat.__get(), CASM::CART);
    casm_coord.voronoi_within();
    return casm_coord.cart();
}

CoordinateEquals_f::CoordinateEquals_f(double tol) : tol(tol) {}
bool CoordinateEquals_f::operator()(const Eigen::Vector3d& ref, const Eigen::Vector3d& other) const
{
    return ref.isApprox(other, tol);
}

} // namespace xtal
} // namespace casmutils
