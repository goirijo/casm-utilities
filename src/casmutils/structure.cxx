#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/lattice.hpp"
#include "casmutils/misc.hpp"
#include <boost/filesystem.hpp>

namespace Extend
{
} // namespace Extend

namespace Rewrap
{
Structure::Structure(const CasmStructure& init_struc) : CasmStructure(init_struc) {}
Structure::Structure(Rewrap::fs::path& filename) : CasmStructure(filename) {}

Structure::Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis)
    : CasmStructure(init_lat)
{
    throw UtilExcept::NotImplemented();
    /* for (const auto& site : init_basis) */
    /* { */
    /*     this->basis.push_back(site.__get()); */
    /*     basis.back().set_lattice(this->lattice(), CASM::CART); */
    /* } */
}

Structure Structure::from_poscar(const fs::path& poscar_path) { return Rewrap::Structure(CasmStructure(poscar_path)); }

bool Structure::is_primitive() const { 
    throw UtilExcept::NotImplemented();
    //return CasmStructure::is_primitive(); 
}

Lattice Structure::lattice() const { return Lattice(this->__get().lattice()); }

void Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode)
{
    CasmStructure* base = this;
    base->set_lattice(new_lattice.__get(), mode);
    return;
}

std::vector<Site> Structure::basis_sites() const
{
    throw UtilExcept::NotImplemented();
    // You are dealing with CASM::Array, but we want to return a std::vector
    //std::vector<Site> basis(this->basis.begin(), this->basis.end());
    //return basis;
}

//**************************************************************************************************************//

Site::operator Coordinate() const { return Coordinate(this->casm_site); }

Site::Site(const Rewrap::Coordinate& init_coord, const std::vector<std::string>& allowed_occupants)
    : casm_site(
          Extend::atomic_site(CASM::xtal::Coordinate(init_coord.cart(), CASM::xtal::Lattice(), CASM::CART), allowed_occupants))
{
    throw UtilExcept::NotImplemented();

    // Avoid an unitialized state.
    /* this->casm_site.set_occ_value(0); */
}

Site::Site(const Eigen::Vector3d& init_coord, const std::vector<std::string>& allowed_occupants)
    : Site(Rewrap::Coordinate(init_coord), allowed_occupants)
{
}

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const Rewrap::Lattice& ref_lattice) const
{
    const_cast<CASM::xtal::Site*>(&this->casm_site)->set_lattice(ref_lattice, CASM::CART);
    return this->casm_site.frac();
}

//**************************************************************************************************************//

Coordinate Coordinate::from_fractional(const Eigen::Vector3d& frac_coord, const Rewrap::Lattice& lat)
{
    CASM::xtal::Coordinate coord(frac_coord, lat, CASM::FRAC);
    return Coordinate(coord);
}

Coordinate Coordinate::from_fractional(double x, double y, double z, const Rewrap::Lattice& lat)
{
    return Coordinate::from_fractional(Eigen::Vector3d(x, y, z), lat);
}

void Coordinate::bring_within(const Lattice& lat)
{
    this->casm_coord.set_lattice(lat, CASM::CART);
    this->casm_coord.within();
    return;
}

Eigen::Vector3d Coordinate::cart() const { return this->casm_coord.cart(); }

Eigen::Vector3d Coordinate::frac(const Rewrap::Lattice& ref_lattice) const
{
    const_cast<CASM::xtal::Coordinate*>(&this->casm_coord)->set_lattice(ref_lattice, CASM::CART);
    return this->casm_coord.frac();
}

Coordinate& Coordinate::operator+=(const Coordinate& coord_to_add)
{
    this->casm_coord += coord_to_add.casm_coord;
    return *this;
}

bool Coordinate::operator==(const Coordinate& coord_to_compare)
{
    return this->casm_coord == coord_to_compare.casm_coord;
}

} // namespace Rewrap
