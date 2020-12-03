#include "../../../autotools.hh"
#include "casmutils/definitions.hpp"
#include <algorithm>
#include <casmutils/misc.hpp>
#include <casmutils/mush/slab.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <tuple>

namespace cu = casmutils;

class SlicingTest : public testing::Test
{
protected:
    using Structure = cu::xtal::Structure;
    void SetUp() override
    {
        hcp_ptr.reset(new Structure(Structure::from_poscar(cu::autotools::input_filesdir / "hcp.vasp")));
        fcc_ptr.reset(new Structure(Structure::from_poscar(cu::autotools::input_filesdir / "fcc.vasp")));
        b2_ptr.reset(new Structure(Structure::from_poscar(cu::autotools::input_filesdir / "b2.vasp")));

        millers_001 << 0, 0, 1;
        millers_010 << 0, 1, 0;
        millers_101 << 1, 0, 1;
        millers_211 << 2, 1, 1;
        millers_111 << 1, 1, 1;
        millers_11bar1 << 1, -1, 1;
    }

    auto make_sliced_structures(const Eigen::Vector3i& miller_indexes)
    {
        return std::make_tuple(cu::xtal::slice_along_plane(*hcp_ptr, miller_indexes),
                               cu::xtal::slice_along_plane(*fcc_ptr, miller_indexes),
                               cu::xtal::slice_along_plane(*b2_ptr, miller_indexes));
    }

    auto make_sliced_lattices(const Eigen::Vector3i& miller_indexes)
    {
        return std::make_tuple(cu::xtal::slice_along_plane(hcp_ptr->lattice(), miller_indexes),
                               cu::xtal::slice_along_plane(fcc_ptr->lattice(), miller_indexes),
                               cu::xtal::slice_along_plane(b2_ptr->lattice(), miller_indexes));
    }

    std::unique_ptr<Structure> hcp_ptr;
    std::unique_ptr<Structure> fcc_ptr;
    std::unique_ptr<Structure> b2_ptr;

    Eigen::Vector3i millers_001;
    Eigen::Vector3i millers_010;
    Eigen::Vector3i millers_101;
    Eigen::Vector3i millers_211;
    Eigen::Vector3i millers_111;
    Eigen::Vector3i millers_11bar1;

    double tol = 1e-5;
};

TEST_F(SlicingTest, SliceLatticeConsistency)
{
    auto miller_set = {millers_001, millers_010, millers_101, millers_211, millers_111, millers_11bar1};
    for (const Eigen::Vector3i& millers : miller_set)
    {
        auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers);
        auto [hcp_slice, fcc_slice, b2_slice] = make_sliced_structures(millers);

        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_slice.lattice(), tol));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_slice.lattice(), tol));
        EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_slice.lattice(), tol));
    }
}

TEST_F(SlicingTest, LatticeSlice001)
{
    auto [hcp_lat, fcc_lat, b2_lat] = make_sliced_lattices(millers_001);
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(hcp_lat, hcp_ptr->lattice(), tol));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_lat, fcc_ptr->lattice(), tol));
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(b2_lat, b2_ptr->lattice(), tol));
}

TEST_F(SlicingTest, LatticeSlice_fcc_111)
{
    cu::xtal::Lattice fcc_slice_lat = cu::xtal::slice_along_plane(fcc_ptr->lattice(), millers_111);
    Eigen::Matrix3i prim_to_3layer_fcc_mat;
    prim_to_3layer_fcc_mat << -1, 0, 1, 1, -1, 1, 0, 1, 1;

    Structure sliced_struc = cu::xtal::slice_along_plane(*fcc_ptr, millers_111);

    cu::xtal::Lattice fcc_3layer = cu::xtal::make_superlattice(fcc_ptr->lattice(), prim_to_3layer_fcc_mat);
    EXPECT_TRUE(cu::is_equal<cu::xtal::LatticeEquals_f>(fcc_slice_lat, fcc_3layer, tol));
}

TEST_F(SlicingTest, StructureSlice_b2_101)
{
    Structure b2_slice = cu::xtal::slice_along_plane(*b2_ptr, millers_101);

    EXPECT_EQ(b2_slice.basis_sites().size(), 4);

    cu::xtal::Site site0(cu::xtal::fractional_to_cartesian(Eigen::Vector3d(0, 0, 0), b2_slice.lattice()), "A");
    cu::xtal::Site site1(cu::xtal::fractional_to_cartesian(Eigen::Vector3d(0.5, 0.5, 0), b2_slice.lattice()), "B");

    cu::xtal::SiteEquals_f is_equal_site(tol);
    EXPECT_TRUE(
        std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), [&](const cu::xtal::Site& site) {
            return is_equal_site(site, site0);
        }) != b2_slice.basis_sites().end());

    EXPECT_TRUE(
        std::find_if(b2_slice.basis_sites().begin(), b2_slice.basis_sites().end(), [&](const cu::xtal::Site& site) {
            return is_equal_site(site, site1);
        }) != b2_slice.basis_sites().end());
}

//******************************************************************************//

class SlabTest : public testing::Test
{
protected:
    using Structure = cu::xtal::Structure;
    void SetUp() override
    {
        b2_101_ptr.reset(new Structure(Structure::from_poscar(cu::autotools::input_filesdir / "b2_101.vasp")));
        b2_101_stack5_ptr.reset(new Structure(cu::mush::make_stacked_slab(*b2_101_ptr, slab_size)));

        cu::xtal::print_poscar(*b2_101_ptr, std::cout);
        cu::xtal::print_poscar(*b2_101_stack5_ptr, std::cout);
    }

    int slab_size = 5;
    std::unique_ptr<Structure> b2_101_ptr;
    std::unique_ptr<Structure> b2_101_stack5_ptr;
};

TEST_F(SlabTest, Consistent_AB_Vectors)
{

    EXPECT_EQ(b2_101_ptr->lattice().a(), b2_101_stack5_ptr->lattice().a());
    EXPECT_EQ(b2_101_ptr->lattice().b(), b2_101_stack5_ptr->lattice().b());
}

TEST_F(SlabTest, FlooredStructure)
{
    auto coord_is_origin = [](const cu::xtal::Site& s) {
        return s.cart().isZero();
    };
    const auto& b2_101_basis = this->b2_101_ptr->basis_sites();
    int ix = 0;
    for (; ix < b2_101_basis.size(); ++ix)
    {
        if (coord_is_origin(b2_101_basis[ix]))
        {
            break;
        }
    }

    std::string original_origin_label = b2_101_basis[ix].label();
    int floor_ix = 0;
    for (; floor_ix < b2_101_basis.size(); ++floor_ix)
    {
        if (b2_101_basis[floor_ix].label() != original_origin_label)
        {
            break;
        }
    }

    Structure floored_structure = cu::mush::make_floored_structure(*(this->b2_101_ptr), floor_ix);
    cu::xtal::print_poscar(floored_structure, std::cout);
    const auto& floored_basis = floored_structure.basis_sites();
    auto origin_site_it = std::find_if(floored_basis.begin(), floored_basis.end(), coord_is_origin);
    EXPECT_TRUE(origin_site_it->label() != original_origin_label);
}

TEST_F(SlabTest, Consistent_C_Vector)
{
    EXPECT_EQ(b2_101_ptr->lattice().c() * slab_size, b2_101_stack5_ptr->lattice().c());
}

TEST_F(SlabTest, OrthogonalSlab)
{
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, slab_size;
    auto simple_stack = cu::xtal::make_superstructure(*b2_101_ptr, stack_mat);

    cu::xtal::LatticeIsEquivalent_f equivalent(1e-8);
    EXPECT_TRUE(equivalent(b2_101_stack5_ptr->lattice(), simple_stack.lattice()));
    EXPECT_TRUE(cu::almost_equal(b2_101_stack5_ptr->lattice().a(), simple_stack.lattice().a(), 1e-8));
    EXPECT_TRUE(cu::almost_equal(b2_101_stack5_ptr->lattice().b(), simple_stack.lattice().b(), 1e-8));
    // Only because c is orthogonal! otherwise it'd probably be false
    EXPECT_TRUE(cu::almost_equal(b2_101_stack5_ptr->lattice().c(), simple_stack.lattice().c(), 1e-8));
}

class OrthogonalizeCVector : public testing::Test
{
protected:
    using Lattice = cu::xtal::Lattice;
    void SetUp() override
    {
        Eigen::Matrix3d col_lat_mat;
        col_lat_mat << 4, 0, 1, 0, 4, 1, 0, 0, 3.2;

        Eigen::AngleAxisd roll(2.5, Eigen::Vector3d::UnitZ());
        Eigen::AngleAxisd yaw(61.9, Eigen::Vector3d::UnitY());
        Eigen::AngleAxisd pitch(-8.3, Eigen::Vector3d::UnitX());
        Eigen::Quaternion<double> q = roll * yaw * pitch;
        rotation_mat = q.matrix();

        col_lat_mat = rotation_mat * col_lat_mat;
        quarter_slant_ptr.reset(new Lattice(col_lat_mat.col(0), col_lat_mat.col(1), col_lat_mat.col(2)));
    }

    void EXPECT_xy_values_for_c_vector(const Eigen::Vector3d& c, double x, double y)
    {
        EXPECT_TRUE(cu::almost_equal(c(0), x, 1e-8));
        EXPECT_TRUE(cu::almost_equal(c(1), y, 1e-8));
    }

    // c vector  slants to 0.25,0.25 fractional coordinates
    std::unique_ptr<Lattice> quarter_slant_ptr;

    // rotate an arbitrary direction so that the vectors aren't aligned in
    // an "easy" direction
    Eigen::Matrix3d rotation_mat;
};

TEST_F(OrthogonalizeCVector, BackToOrigin)
{
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, 4;

    auto stack = cu::xtal::make_superlattice(*quarter_slant_ptr, stack_mat);
    auto ortho = cu::mush::make_aligned(cu::mush::orthogonalize_c_vector(stack));
    EXPECT_xy_values_for_c_vector(ortho.c(), 0, 0);
}

TEST_F(OrthogonalizeCVector, BackwardQuarter)
{
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, 3;

    auto stack = cu::xtal::make_superlattice(*quarter_slant_ptr, stack_mat);
    auto ortho = cu::mush::make_aligned(cu::mush::orthogonalize_c_vector(stack));
    EXPECT_xy_values_for_c_vector(ortho.c(), -1, -1);
}

TEST_F(OrthogonalizeCVector, ForwardQuarter)
{
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, 9;

    auto stack = cu::xtal::make_superlattice(*quarter_slant_ptr, stack_mat);
    auto ortho = cu::mush::make_aligned(cu::mush::orthogonalize_c_vector(stack));
    EXPECT_xy_values_for_c_vector(ortho.c(), 1, 1);
}

TEST_F(OrthogonalizeCVector, EquivalentToStack)
{
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, 7;

    auto stack = cu::xtal::make_superlattice(*quarter_slant_ptr, stack_mat);
    auto ortho = cu::mush::orthogonalize_c_vector(stack);

    cu::xtal::LatticeIsEquivalent_f equivalent(1e-8);
    EXPECT_TRUE(equivalent(stack, ortho));

    // TODO: make this a binary comparator ffs
    cu::xtal::LatticeEquals_f equal(1e-8);
    EXPECT_FALSE(equal(stack, ortho));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
