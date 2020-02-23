#include <casmutils/xtal/structure.hpp>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/definitions.hpp>

namespace Extend
{
} // namespace Extend

namespace Rewrap
{
Structure::Structure(const CASM::xtal::BasicStructure& init_struc) : structure_lattice(init_struc.lattice())
{
    throw UtilExcept::NotImplemented();
}

Structure::Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis)
    : structure_lattice(init_lat)
{
    throw UtilExcept::NotImplemented();
    /* for (const auto& site : init_basis) */
    /* { */
    /*     this->basis.push_back(site.__get()); */
    /*     basis.back().set_lattice(this->lattice(), CASM::CART); */
    /* } */
}

Structure Structure::from_poscar(const fs::path& poscar_path) { throw UtilExcept::NotImplemented(); }

const Lattice& Structure::lattice() const { return this->structure_lattice; }

void Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode)
{
    throw UtilExcept::NotImplemented();
    return;
}

const std::vector<Site>& Structure::basis_sites() const
{
    throw UtilExcept::NotImplemented();
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

} // namespace Rewrap
