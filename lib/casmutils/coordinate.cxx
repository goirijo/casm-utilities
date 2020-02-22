#include "casmutils/coordinate.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/lattice.hpp"
#include "casmutils/misc.hpp"

namespace Rewrap
{

Coordinate Coordinate::from_fractional(const Eigen::Vector3d& frac_coord, const Rewrap::Lattice& lat)
{
    CASM::xtal::Coordinate coord(frac_coord, lat.__get(), CASM::FRAC);
    return Coordinate(coord);
}

Coordinate Coordinate::from_fractional(double x, double y, double z, const Rewrap::Lattice& lat)
{
    return Coordinate::from_fractional(Eigen::Vector3d(x, y, z), lat);
}

void Coordinate::bring_within(const Lattice& lat)
{
    this->casm_coord.set_lattice(lat.__get(), CASM::CART);
    this->casm_coord.within();
    return;
}

Eigen::Vector3d Coordinate::cart() const { return this->casm_coord.cart(); }

Eigen::Vector3d Coordinate::frac(const Rewrap::Lattice& ref_lattice) const
{
    const_cast<CASM::xtal::Coordinate*>(&this->casm_coord)->set_lattice(ref_lattice.__get(), CASM::CART);
    return this->casm_coord.frac();
}

Coordinate& Coordinate::operator+=(const Coordinate& coord_to_add)
{
    this->casm_coord += coord_to_add.casm_coord;
    return *this;
}

bool Coordinate::operator==(const Coordinate& coord_to_compare)
{
    return this->casm_coord == coord_to_compare.casm_coord;
}

} // namespace Rewrap
