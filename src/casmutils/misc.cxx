#include <casm/clex/PrimClex.hh>

namespace Extend
{
CASM::PrimClex quiet_primclex(CASM::Structure& prim)
{
    CASM::Log log(std::cout, 0);
    CASM::Logging logging(log);
    CASM::PrimClex pclex(prim, logging);
    return pclex;
}

CASM::Site atomic_site(const CASM::Coordinate& coord, const std::vector<std::string>& allowed_species)
{
    CASM::Array<CASM::Molecule> allowed_molecules;
    for (auto specie : allowed_species)
    {
        // make_atom returns a molecule lol
        allowed_molecules.push_back(make_atom(specie, coord.home()));
    }

    CASM::Site site(coord, "");
    site.set_site_occupant(allowed_molecules);

    return site;
}
} // namespace Extend

namespace IO
{
Eigen::IOFormat coord_format() { return Eigen::IOFormat(4, 0, ", ", "\n", "[", "]"); }
} // namespace IO
