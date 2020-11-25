#ifndef UTILS_SITE_HH
#define UTILS_SITE_HH

#include <casm/crystallography/Site.hh>
#include <casmutils/definitions.hpp>

namespace casmutils
{
namespace xtal
{

class Lattice;
class Coordinate;

/**
 * A coordinate and type of species. Even though it's implemented
 * with CASM::xtal::Site, it can *NOT* describe mutiple possible occupants.
 */

class Site
{
public:
    Site() = delete;
    Site(const CASM::xtal::Site& init_site, int occupant);
    Site(const Coordinate& init_coord, const std::string& occupant_name);
    Site(const Eigen::Vector3d& init_coord, const std::string& occupant_name);

    /// Allow casting to Coordinate, by stripping everything away except the Cartesian position
    operator Coordinate() const;

    /// Retreive the Cartesian values of the coordinate
    Eigen::Vector3d cart() const;
    /// Retreive the fractional values of the coordinate relative to the provided lattice
    Eigen::Vector3d frac(const Lattice& ref_lattice) const;

    /// Name of the species residing on the site
    std::string label() const;

    /// Access  the CASM implementation within.
    /* const CASM::xtal::Site& __get() const {return casm_site;}; */

private:
    CASM::xtal::Site casm_site;
};

/// This functor class provides a unary equals operator
/// for casmutils::xtal::Site
class SiteEquals_f
{
public:
    /// Determines whether test is equal to reference site with a tolerance
    SiteEquals_f(double tol);
    /// Returns true if other is equal to ref_site
    bool operator()(const Site& ref_site, const Site& other) const;

private:
    double tol;
};
} // namespace xtal
} // namespace casmutils

#endif
