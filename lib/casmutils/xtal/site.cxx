#include "casmutils/xtal/site.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/lattice.hpp"
#include <boost/filesystem.hpp>

namespace extend
{
} // namespace extend

namespace rewrap
{

Site::operator Coordinate() const { return Coordinate(this->casm_site); }

Site::Site(const rewrap::Coordinate& init_coord, const std::string& occupant_name)
    : casm_site(CASM::xtal::Site(init_coord.__get(), occupant_name))

{
    // Avoid a multioccupant site
    assert(casm_site.allowed_occupants().size() == 1);
    // Avoid an unitialized state.
    this->casm_site.set_label(0);
}

std::string Site::label() const { return this->casm_site.allowed_occupants()[this->casm_site.label()]; }

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const rewrap::Lattice& ref_lattice) const
{
    return ref_lattice.column_vector_matrix().inverse() * this->cart();
}

} // namespace rewrap
