#include "casmutils/xtal/site.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/lattice.hpp"
#include <boost/filesystem.hpp>

namespace extend
{
} // namespace extend

namespace casmutils
{
namespace xtal
{

Site::operator Coordinate() const { return Coordinate(this->casm_site); }

Site::Site(const Coordinate& init_coord, const std::string& occupant_name)
    : casm_site(CASM::xtal::Site(init_coord.__get(), occupant_name))

{
    // Avoid a multioccupant site
    assert(casm_site.allowed_occupants().size() == 1);
    // Avoid an unitialized state.
    this->casm_site.set_label(0);
}
Site::Site(const CASM::xtal::Site& init_site, int occupant) : casm_site(init_site)
{
    // Avoid an unitialized state.
    this->casm_site.set_label(occupant);
}
std::string Site::label() const { return this->casm_site.allowed_occupants()[this->casm_site.label()]; }

Eigen::Vector3d Site::cart() const { return this->casm_site.cart(); }

Eigen::Vector3d Site::frac(const Lattice& ref_lattice) const
{
    return ref_lattice.column_vector_matrix().inverse() * this->cart();
}

SiteEquals_f::SiteEquals_f(const Site& ref_site, double tol) : ref_site(ref_site), tol(tol) {}
bool SiteEquals_f::operator()(const Site& other)
{
    return is_equal<CoordinateEquals_f>(static_cast<Coordinate>(ref_site), static_cast<Coordinate>(other), tol) &&
           ref_site.label() == other.label();
}

} // namespace xtal
} // namespace casmutils
