#ifndef UTILS_COORDINATE_HH
#define UTILS_COORDINATE_HH

#include <casm/crystallography/Coordinate.hh>
#include <casmutils/definitions.hpp>

namespace casmutils
{
namespace xtal
{

class Lattice;

// TODO: Just erase this class, it's super annoying. You can do all of this
// stuff with functions

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
    static Coordinate from_fractional(const Eigen::Vector3d& frac_coord, const Lattice& lat);
    static Coordinate from_fractional(double x, double y, double z, const Lattice& lat);

    /// Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    /// Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const Lattice& ref_lattice) const;

    Coordinate operator+(const Coordinate& coord_to_add) const;
    /// Translate *this by the given coordinate
    Coordinate& operator+=(const Coordinate& coord_to_add);

    // TODO: These are annoying as member functions. Make them standalone.
    // and let them accept Eigen::Vector3d as well. Same with the wigner-seitz
    /// Bring *this within the given lattice
    void bring_within(const Lattice& lat);

    /// Bring *this within the given lattice
    [[nodiscard]] Coordinate bring_within(const Lattice& lat) const;

    /// Bring *this within the given lattice
    void bring_within_wigner_seitz(const Lattice& lat);

    /// Bring *this within the Wigner-Seitz cell of the given lattice
    [[nodiscard]] Coordinate bring_within_wigner_seitz(const Lattice& lat) const;

    /// Access  the CASM implementation within.
    const CASM::xtal::Coordinate& __get() const { return casm_coord; };

private:
    /// Use the CASM implementation to forward any functionality you want
    CASM::xtal::Coordinate casm_coord;
};

// TODO: make this binary, not unary (same with all other comparators)
struct CoordinateEquals_f
{
    /// for casmutils::xtal::Coordinate
public:
    /// Determines whether test is equal to reference site with a tolerance
    CoordinateEquals_f(const Coordinate& ref_coordinate, double tol);
    /// Returns true if other is equal to the Coordinate the comparator was constructed with
    bool operator()(const Coordinate& other);

private:
    Coordinate ref_coordinate;
    double tol;
};

} // namespace xtal
} // namespace casmutils

#endif
