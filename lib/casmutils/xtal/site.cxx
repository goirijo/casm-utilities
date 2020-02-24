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
    throw except::NotImplemented();

    // Avoid an unitialized state.
    /* this->casm_site.set_occ_value(0); */
}

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const rewrap::Lattice& ref_lattice) const { throw except::NotImplemented(); }

} // namespace rewrap
