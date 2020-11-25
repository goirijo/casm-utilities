#include "casmutils/xtal/site.hpp"
#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/lattice.hpp"

namespace extend
{
} // namespace extend

namespace casmutils
{
namespace xtal
{

Site::Site(const Eigen::Vector3d& init_coord, const std::string& occupant_name)
    : casm_site(CASM::xtal::Site(CASM::xtal::Coordinate(init_coord, CASM::xtal::Lattice(), CASM::CART), occupant_name))
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

SiteEquals_f::SiteEquals_f(double tol) : tol(tol) {}
bool SiteEquals_f::operator()(const Site& ref_site, const Site& other) const
{
    return is_equal<CoordinateEquals_f>(ref_site.cart(), other.cart(), tol) && ref_site.label() == other.label();
}

} // namespace xtal
} // namespace casmutils
