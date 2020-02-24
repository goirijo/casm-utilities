#include "casmutils/misc.hpp"
#include "casmutils/exceptions.hpp"
#include <casm/clex/PrimClex.hh>
#include <casm/crystallography/Site.hh>
#include <casm/external/Eigen/Core>

namespace extend
{
CASM::PrimClex quiet_primclex(CASM::Structure& prim)
{
    throw except::NotImplemented();
    /* CASM::Log log(std::cout, 0); */
    /* CASM::Logging logging(log); */
    /* return CASM::PrimClex(prim, logging); */
}

CASM::xtal::Site atomic_site(const CASM::xtal::Coordinate& coord, const std::vector<std::string>& allowed_species)
{
    throw except::NotImplemented(); // PS: This is already a static method of Molecule
    /* CASM::Array<CASM::Molecule> allowed_molecules; */
    /* for (auto specie : allowed_species) */
    /* { */
    /*     // make_atom returns a molecule lol */
    /*     allowed_molecules.push_back(make_atom(specie, coord.home())); */
    /* } */

    /* CASM::Site site(coord, ""); */
    /* site.set_site_occupant(allowed_molecules); */

    /* return site; */
}
} // namespace extend

namespace io
{
Eigen::IOFormat coord_format() { return Eigen::IOFormat(4, 0, 0, ", ", "\n", "[", "]"); }
} // namespace io
