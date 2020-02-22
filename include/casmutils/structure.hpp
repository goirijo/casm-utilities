#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casmutils/lattice.hpp"
#include "definitions.hpp"
#include <iostream>

namespace Rewrap
{
class Lattice;

/**
 * The Rewrap version of Coordinate does *not* have a home Lattice,
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
    static Coordinate from_fractional(const Eigen::Vector3d& frac_coord, const Rewrap::Lattice& lat);
    static Coordinate from_fractional(double x, double y, double z, const Rewrap::Lattice& lat);

    /// Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    /// Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const Rewrap::Lattice& ref_lattice) const;

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

/**
 * A coordinate and type of species. Even though it's implemented
 * with CASM::xtal::Site, it can *NOT* describe mutiple possible occupants.
 */

class Site
{
public:
    Site() = delete;
    Site(const CASM::xtal::Site& init_site, int occupant) : casm_site(init_site) {}
    Site(const Coordinate& init_coord, const std::string& occupant_name);
    /* Site(const Eigen::Vector3d& init_coord, const std::string& occupant_name); */

    /// Allow casting to Coordinate, by stripping everything away except the Cartesian position
    operator Coordinate() const;

    /// Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    /// Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const Rewrap::Lattice& ref_lattice) const;

    /// Name of the species residing on the site
    std::string label() const;

    /// Access  the CASM implementation within.
    /* const CASM::xtal::Site& __get() const {return casm_site;}; */

private:
    CASM::xtal::Site casm_site;
};

/**
 * Describes a current, fixed state of a crystal. Composed of a lattice,
 * and collection of basis atoms (Sites). Each site contains the position of
 * the species, and a label for the species.
 */

class Structure
{
public:
    Structure() = delete;

    /// Construct by providing a path to a POSCAR like file
    static Structure from_poscar(const fs::path& poscar_path);

    /// Construct from parent CASM class
    Structure(const CASM::xtal::SimpleStructure& init_struc);
    Structure(const CASM::xtal::BasicStructure& init_struc);

    /// Construct with a lattice and list of sites (basis)
    Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis);

    /// Returns a copy of the current lattice of the structure
    const Lattice& lattice() const;

    /// Give the structure a new lattice, and either keep the Cartesian, or fractional coordinates of the basis
    void set_lattice(const Lattice& new_lattice, COORD_TYPE mode);

    /// Give the structure a new lattice, and either keep the Cartesian, or fractional coordinates of the basis
    Structure set_lattice(const Lattice& new_lattice, COORD_TYPE mode) const;

    /// Return a copy of all the basis sites
    const std::vector<Site>& basis_sites() const;

    /// Retreive the CASM implementations of *this
    template <typename CASMType> const CASMType& __get() const;

private:
    /// CASM::SimpleStructure representation
    CASM::xtal::SimpleStructure casm_simplestructure;

    /// CASM::SimpleStructure representation
    CASM::xtal::BasicStructure casm_basicstructure;

    /// Rewrap representation of the basis
    std::vector<Site> basis;

    /// Rewrap representation of the lattice
    Lattice structure_lattice;
};
} // namespace Rewrap

#endif
