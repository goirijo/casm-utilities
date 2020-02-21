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
Structure::Structure(const CASM::xtal::BasicStructure& init_struc): structure_lattice(init_struc.lattice()){
    throw UtilExcept::NotImplemented();
}

Structure::Structure(const Rewrap::Lattice& init_lat, const std::vector<Rewrap::Site>& init_basis):structure_lattice(init_lat)
{
    throw UtilExcept::NotImplemented();
    /* for (const auto& site : init_basis) */
    /* { */
    /*     this->basis.push_back(site.__get()); */
    /*     basis.back().set_lattice(this->lattice(), CASM::CART); */
    /* } */
}

Structure Structure::from_poscar(const fs::path &poscar_path){
    throw UtilExcept::NotImplemented();
}

const Lattice& Structure::lattice() const { return this->structure_lattice;}

void Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode)
{
    throw UtilExcept::NotImplemented();
    return;
}

const std::vector<Site>& Structure::basis_sites() const
{
    throw UtilExcept::NotImplemented();
    // You are dealing with CASM::Array, but we want to return a std::vector
    //std::vector<Site> basis(this->basis.begin(), this->basis.end());
    //return basis;
}

    ///Return *this as a CASM::BasicStructure
    template<>
    const CASM::xtal::SimpleStructure& Structure::__get<CASM::xtal::SimpleStructure>() const {return this->casm_simplestructure;}

    ///Return *this as a CASM::BasicStructure
    template<>
    const CASM::xtal::BasicStructure& Structure::__get<CASM::xtal::BasicStructure>() const {return this->casm_basicstructure;}


//**************************************************************************************************************//

Site::operator Coordinate() const { return Coordinate(this->casm_site); }

Site::Site(const Rewrap::Coordinate& init_coord, const std::string& occupant_name)
    : casm_site(CASM::xtal::Site(init_coord.__get(),occupant_name))
          
{
    throw UtilExcept::NotImplemented();

    // Avoid an unitialized state.
    /* this->casm_site.set_occ_value(0); */
}

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const Rewrap::Lattice& ref_lattice) const
{
    throw UtilExcept::NotImplemented();
}

//**************************************************************************************************************//

Coordinate Coordinate::from_fractional(const Eigen::Vector3d& frac_coord, const Rewrap::Lattice& lat)
{
    CASM::xtal::Coordinate coord(frac_coord, lat.__get(), CASM::FRAC);
    return Coordinate(coord);
}

Coordinate Coordinate::from_fractional(double x, double y, double z, const Rewrap::Lattice& lat)
{
    return Coordinate::from_fractional(Eigen::Vector3d(x, y, z), lat);
}

void Coordinate::bring_within(const Lattice& lat)
{
    this->casm_coord.set_lattice(lat.__get(), CASM::CART);
    this->casm_coord.within();
    return;
}

Eigen::Vector3d Coordinate::cart() const { return this->casm_coord.cart(); }

Eigen::Vector3d Coordinate::frac(const Rewrap::Lattice& ref_lattice) const
{
    const_cast<CASM::xtal::Coordinate*>(&this->casm_coord)->set_lattice(ref_lattice.__get(), CASM::CART);
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
