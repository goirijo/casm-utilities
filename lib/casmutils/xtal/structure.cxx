#include <casmutils/definitions.hpp>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/structure.hpp>

namespace extend
{
} // namespace extend

namespace rewrap
{
Structure::Structure(const CASM::xtal::BasicStructure& init_struc) : structure_lattice(init_struc.lattice())
{
    throw except::NotImplemented();
}

Structure::Structure(const rewrap::Lattice& init_lat, const std::vector<rewrap::Site>& init_basis)
    : structure_lattice(init_lat)
{
    throw except::NotImplemented();
    /* for (const auto& site : init_basis) */
    /* { */
    /*     this->basis.push_back(site.__get()); */
    /*     basis.back().set_lattice(this->lattice(), CASM::CART); */
    /* } */
}

Structure Structure::from_poscar(const fs::path& poscar_path) { throw except::NotImplemented(); }

const Lattice& Structure::lattice() const { return this->structure_lattice; }

void Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode)
{
    throw except::NotImplemented();
    return;
}

[[nodiscard]] Structure Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode) const{
	throw except::NotImplemented();
}
const std::vector<Site>& Structure::basis_sites() const
{
    throw except::NotImplemented();
    // You are dealing with CASM::Array, but we want to return a std::vector
    // std::vector<Site> basis(this->basis.begin(), this->basis.end());
    // return basis;
}



/// Return *this as a CASM::BasicStructure
template <> const CASM::xtal::SimpleStructure& Structure::__get<CASM::xtal::SimpleStructure>() const
{
    return this->casm_simplestructure;
}

/// Return *this as a CASM::BasicStructure
template <> const CASM::xtal::BasicStructure& Structure::__get<CASM::xtal::BasicStructure>() const
{
    return this->casm_basicstructure;
}

} // namespace rewrap
