#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
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
    Coordinate(const CASM::Coordinate& init_coord) : casm_coord(init_coord) {}
    Coordinate(const Eigen::Vector3d& cart_coord)
        : Coordinate(CASM::Coordinate(cart_coord, CASM::Lattice(), CASM::CART)){};
    Coordinate(double x, double y, double z)
        : Coordinate(Eigen::Vector3d(x,y,z)){}

    ///Initialize the Coordinate by providing fractional coordinates and the corresponding lattice
    static Coordinate from_fractional(const Eigen::Vector3d& frac_coord, const Rewrap::Lattice& lat);
    static Coordinate from_fractional(double x, double y, double z, const Rewrap::Lattice& lat);

    ///Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    ///Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const Rewrap::Lattice& ref_lattice) const;

    ///Translate *this by the given coordinate
    Coordinate &operator+=(const Coordinate &coord_to_add);

    ///Returns true if the distance between the coordinates is within CASM::TOL
    bool operator==(const Coordinate &coord_to_compare);

private:
    /// Use the CASM implementation to forward any functionality you want
    CASM::Coordinate casm_coord;
};

class Site
{
public:
    Site() = delete;
    Site(const CASM::Site& init_site): casm_site(init_site){}
    Site(Coordinate& init_coord, const std::vector<std::string>& allowed_occupants);

    ///Allow casting to Coordinate, by stripping everything away except the Cartesian position
    operator Coordinate() const;

    ///Access  the CASM implementation within.
    const CASM::Site& __get() const {return casm_site;};

private:
    CASM::Site casm_site;
};

class Structure : public CASM::Structure
{
public:
    Structure() = delete;

    /// Construct by providing a path to a POSCAR like file
    static Structure from_poscar(const fs::path& poscar_path);

    /// Construct from parent CASM class
    Structure(const CASM::Structure& init_struc);
    Structure(Rewrap::fs::path& filename);

    /// Construct with a lattice and list of sites (basis)
    Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis);

    /// Returns true if the structure is already primitive
    bool is_primitive() const;

    /// Add a new site to the basis
    /* void push_basis(const Rewrap::Site& new_basis_site); */

    ///Return a copy of all the basis sites
    std::vector<Site> basis_sites() const;

private:
};
} // namespace Rewrap


#endif
