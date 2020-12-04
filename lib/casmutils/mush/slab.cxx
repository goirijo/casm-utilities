#include "casmutils/xtal/structure.hpp"
#include <casmutils/mush/slab.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/frankenstein.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure_tools.hpp>

// Extract to casm-utilities
#include <casm/crystallography/Lattice.hh>
#include <casm/crystallography/Superlattice.hh>
#include <vector>

namespace casmutils
{
namespace mush
{
Eigen::Matrix3d slab_unit_vectors(const xtal::Lattice& slab)
{
    if (slab.column_vector_matrix().determinant() < 0)
    {
        throw std::runtime_error("Encountered a left handed lattice when making slab unit vectors.");
    }

    Eigen::Matrix3d slab_span;
    slab_span.col(0) = slab.a().normalized();
    slab_span.col(2) = slab.a().cross(slab.b()).normalized();
    slab_span.col(1) = slab_span.col(2).cross(slab_span.col(0)).normalized();

    assert(slab_span.determinant() * slab.column_vector_matrix().determinant() > 0);
    assert(almost_equal(std::abs(slab_span.determinant()), 1.0, 1e-10));

    return slab_span;
}

/// Applies the given transformation the the *column* vector matrix representation
/// of the lattice
xtal::Lattice make_transformed_lattice(const xtal::Lattice& lat, const Eigen::Matrix3d& transform)
{
    return xtal::Lattice(transform * lat.column_vector_matrix());
}

Eigen::Matrix3d make_alignment_matrix(const xtal::Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = slab_unit_vectors(lat).inverse();
    return lat_span_to_standard;
}

xtal::Lattice make_aligned(const xtal::Lattice& lat)
{
    Eigen::Matrix3d lat_span_to_standard = make_alignment_matrix(lat);
    xtal::Lattice aligned_lat = make_transformed_lattice(lat, lat_span_to_standard);

    assert(CASM::almost_equal(aligned_lat.a()(2), 0.0, 1e-10));
    assert(CASM::almost_equal(aligned_lat.a()(1), 0.0, 1e-10));
    assert(CASM::almost_equal(aligned_lat.b()(2), 0.0, 1e-10));

    assert(CASM::almost_equal((aligned_lat.a().cross(aligned_lat.b())).normalized(), Eigen::Vector3d(0, 0, 1), 1e-10));
    assert(CASM::almost_equal(
        lat.column_vector_matrix().determinant(), aligned_lat.column_vector_matrix().determinant(), 1e-10));
    return aligned_lat;
}

xtal::Structure make_aligned(xtal::Structure struc)
{
    make_aligned(&struc);
    return struc;
}

void make_aligned(xtal::Structure* struc)
{
    struc->set_lattice(make_aligned(struc->lattice()), xtal::FRAC);
    return;
}

xtal::Lattice orthogonalize_c_vector(const xtal::Lattice& lat)
{
    assert(lat.column_vector_matrix().determinant() > 0);
    Eigen::Vector3d normal = lat.a().cross(lat.b()).normalized();

    Eigen::Vector3d ortho_component = lat.c().dot(normal) * normal;
    Eigen::Vector3d plane_component = lat.c() - ortho_component;

    // In plane component should have no c component

    auto plane_component_frac = xtal::cartesian_to_fractional(plane_component, lat);
    assert(almost_equal(plane_component_frac(2), 0.0, 1e-8));

    plane_component = xtal::bring_within_wigner_seitz(plane_component, lat);
    Eigen::Vector3d final_c = plane_component + ortho_component;
    return xtal::Lattice(lat.a(), lat.b(), final_c);
}

xtal::Structure orthogonalize_c_vector(const xtal::Structure& struc)
{
    return struc.set_lattice(orthogonalize_c_vector(struc.lattice()), xtal::CART);
}

xtal::Structure make_stacked_slab(const xtal::Structure& slab_unit, int stacks)
{
    // Always stack along c-direction
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, stacks;

    xtal::Structure stack = xtal::make_superstructure(slab_unit, stack_mat);
    stack.set_lattice(orthogonalize_c_vector(stack.lattice()), xtal::CART);

    return stack;
}

xtal::Structure make_floored_structure(const xtal::Structure& shiftable_struc, int floor_atom_ix)
{
    // Index 0 means no translation
    if (floor_atom_ix < 1)
    {
        return shiftable_struc;
    }

    const Eigen::Vector3d cart_shift = -shiftable_struc.basis_sites()[floor_atom_ix].cart();
    return frankenstein::translate_basis(shiftable_struc, cart_shift);
}
} // namespace mush
} // namespace casmutils
