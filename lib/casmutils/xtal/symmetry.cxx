#include "casmutils/xtal/structure.hpp"
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/SymTools.hh>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/symmetry.hpp>

namespace casmutils
{
namespace xtal
{
std::vector<sym::CartOp> make_point_group(const Lattice& lat, double tol)
{
    return CASM::xtal::make_point_group(lat.__get(), tol);
}

std::vector<sym::CartOp> make_factor_group(const Structure& struc, double tol)
{
    return CASM::xtal::make_factor_group(struc.__get<CASM::xtal::BasicStructure>(), tol);
}

Lattice symmetrize(const Lattice& noisy_lattice, const std::vector<sym::CartOp>& enforced_point_group)
{
    return CASM::xtal::symmetrize(noisy_lattice.__get(), enforced_point_group);
}

} // namespace xtal
} // namespace casmutils
