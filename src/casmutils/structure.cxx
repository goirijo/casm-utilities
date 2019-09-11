#include <boost/filesystem.hpp>
#include <casmutils/structure.hpp>
#include <casmutils/lattice.hpp>

namespace Rewrap
{
Structure::Structure(const CASM::Structure& init_struc) : CASM::Structure(init_struc) {}
Structure::Structure(Rewrap::fs::path& filename) : CASM::Structure(filename) {}

Structure::Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis):CASM::Structure(init_lat)
{
    for(const auto& site : init_basis)
    {
        this->basis.push_back(site.__get());
    }
}

Structure Structure::from_poscar(const fs::path& poscar_path)
{
    return Rewrap::Structure(CASM::Structure(poscar_path));
}

bool Structure::is_primitive() const { return CASM::Structure::is_primitive(); }

std::vector<Site> Structure::basis_sites() const
{
    //You are dealing with CASM::Array, but we want to return a std::vector
    std::vector<Site> basis(this->basis.begin(), this->basis.end());
    return basis;
}

Site::operator Coordinate() const
{
    return Coordinate(this->casm_site);
}

} // namespace Rewrap
