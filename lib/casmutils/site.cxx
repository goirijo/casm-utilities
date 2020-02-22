#include "casmutils/site.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/lattice.hpp"
#include "casmutils/coordinate.hpp"
#include "casmutils/misc.hpp"
#include <boost/filesystem.hpp>

namespace Extend
{
} // namespace Extend

namespace Rewrap
{

Site::operator Coordinate() const { return Coordinate(this->casm_site); }

Site::Site(const Rewrap::Coordinate& init_coord, const std::string& occupant_name)
    : casm_site(CASM::xtal::Site(init_coord.__get(), occupant_name))

{
    throw UtilExcept::NotImplemented();

    // Avoid an unitialized state.
    /* this->casm_site.set_occ_value(0); */
}

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const Rewrap::Lattice& ref_lattice) const { throw UtilExcept::NotImplemented(); }

} // namespace Rewrap
