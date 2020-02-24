#ifndef UTILS_COORDINATE_HH
#define UTILS_COORDINATE_HH

#include <casm/crystallography/Coordinate.hh>
#include <casmutils/definitions.hpp>

namespace rewrap
{
class Lattice;

/**
 * The rewrap version of Coordinate does *not* have a home Lattice,
 * and instead is always set to CART mode. Any routines that involve
 * a lattice require passing it as an argument.
 */

class Coordinate
{
public:
    Coordinate() = delete;
    Coordinate(const CASM::xtal::Coordinate& init_coord) : casm_coord(init_coord) {}
    Coordinate(const Eigen::Vector3d& cart_coord)
        : Coordinate(CASM::xtal::Coordinate(cart_coord, CASM::xtal::Lattice(), CASM::CART)){};
    Coordinate(double x, double y, double z) : Coordinate(Eigen::Vector3d(x, y, z)) {}

    /// Initialize the Coordinate by providing fractional coordinates and the corresponding lattice
    static Coordinate from_fractional(const Eigen::Vector3d& frac_coord, const rewrap::Lattice& lat);
    static Coordinate from_fractional(double x, double y, double z, const rewrap::Lattice& lat);

    /// Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    /// Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const rewrap::Lattice& ref_lattice) const;

    /// Translate *this by the given coordinate
    Coordinate& operator+=(const Coordinate& coord_to_add);

    /// Bring *this within the given lattice
    void bring_within(const Lattice& lat);

    /// Returns true if the distance between the coordinates is within CASM::TOL
    bool operator==(const Coordinate& coord_to_compare);

    /// Access  the CASM implementation within.
    const CASM::xtal::Coordinate& __get() const { return casm_coord; };

private:
    /// Use the CASM implementation to forward any functionality you want
    CASM::xtal::Coordinate casm_coord;
};
} // namespace rewrap

namespace casmutils
{
namespace xtal
{
using rewrap::Coordinate;
}
} // namespace casmutils

#endif
