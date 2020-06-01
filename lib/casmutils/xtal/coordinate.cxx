#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/lattice.hpp"

namespace casmutils
{
namespace xtal
{
Coordinate Coordinate::from_fractional(const Eigen::Vector3d& frac_coord, const Lattice& lat)
{
    CASM::xtal::Coordinate coord(frac_coord, lat.__get(), CASM::FRAC);
    return Coordinate(coord);
}

Coordinate Coordinate::from_fractional(double x, double y, double z, const Lattice& lat)
{
    return Coordinate::from_fractional(Eigen::Vector3d(x, y, z), lat);
}

void Coordinate::bring_within(const Lattice& lat)
{
    this->casm_coord.set_lattice(lat.__get(), CASM::CART);
    this->casm_coord.within();
    return;
}

Coordinate Coordinate::bring_within(const Lattice& lat) const
{
    Coordinate copy_coord(*this);
    copy_coord.bring_within(lat);
    return copy_coord;
}

void Coordinate::bring_within_wigner_seitz(const Lattice& lat)
{
    this->casm_coord.set_lattice(lat.__get(), CASM::CART);
    this->casm_coord.voronoi_within();
    return;
}

Coordinate Coordinate::bring_within_wigner_seitz(const Lattice& lat) const
{
    Coordinate copy_coord(*this);
    copy_coord.bring_within_wigner_seitz(lat);
    return copy_coord;
}

Eigen::Vector3d Coordinate::cart() const { return this->casm_coord.cart(); }

Eigen::Vector3d Coordinate::frac(const Lattice& ref_lattice) const
{
    const_cast<CASM::xtal::Coordinate*>(&this->casm_coord)->set_lattice(ref_lattice.__get(), CASM::CART);
    return this->casm_coord.frac();
}

Coordinate& Coordinate::operator+=(const Coordinate& coord_to_add)
{
    this->casm_coord += coord_to_add.casm_coord;
    return *this;
}

Coordinate Coordinate::operator+(const Coordinate& coord_to_add) const
{
    Coordinate summed_coord = *this;
    summed_coord += coord_to_add;
    return summed_coord;
}

CoordinateEquals_f::CoordinateEquals_f(const Coordinate& ref_coordinate, double tol)
    : ref_coordinate(ref_coordinate), tol(tol)
{
}
bool CoordinateEquals_f::operator()(const Coordinate& other)
{
    if ((ref_coordinate.cart() - other.cart()).norm() < tol)
    {
        return true;
    }

    return false;
}

} // namespace xtal
} // namespace casmutils
