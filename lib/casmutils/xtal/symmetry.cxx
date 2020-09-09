#include "casmutils/xtal/structure.hpp"
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/SymTools.hh>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
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

Structure symmetrize(const Structure& noisy_structure, const std::vector<sym::CartOp>& enforced_factor_group)
{
    Lattice corrected_lattice = symmetrize(noisy_structure.lattice(), enforced_factor_group);
    Structure structure_with_correct_lattice = noisy_structure;
    structure_with_correct_lattice.set_lattice(corrected_lattice, FRAC);
    return CASM::xtal::symmetrize(structure_with_correct_lattice.__get<CASM::xtal::BasicStructure>(),
                                  enforced_factor_group);
}

Eigen::Vector3d operator*(const sym::CartOp& sym_op, const Eigen::Vector3d& vector3d)
{
    return sym_op.matrix * vector3d + sym_op.translation;
}

Site operator*(const sym::CartOp& sym_op, const Site& site) { return Site{sym_op * site.cart(), site.label()}; }

Coordinate operator*(const sym::CartOp& sym_op, const Coordinate& coordinate)
{
    return Coordinate{sym_op * coordinate.cart()};
}

} // namespace xtal
} // namespace casmutils
