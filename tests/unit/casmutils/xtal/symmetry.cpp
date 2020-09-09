#include "../../../autotools.hh"
#include "casmutils/sym/cartesian.hpp"
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>
#include <vector>

namespace cu = casmutils;

class CrystalGroupTest : public testing::Test
{
protected:
    using Structure = cu::xtal::Structure;

    void SetUp() override
    {
        cu::fs::path primitive_fcc_path(cu::autotools::input_filesdir / "primitive_fcc_Ni.vasp");
        primitive_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(primitive_fcc_path));
    }

    std::unique_ptr<Structure> primitive_fcc_Ni_ptr;
    double tol = 1e-10;
};

class SymmetrizeTest : public testing::Test
{
protected:
    using Lattice = cu::xtal::Lattice;
    using Structure = cu::xtal::Structure;
    using CartOp = cu::sym::CartOp;
    using Coordinate = cu::xtal::Coordinate;
    void SetUp() override
    {
        Eigen::Matrix3d cubic_lat_matrix = 3 * Eigen::Matrix3d::Identity();
        Eigen::Matrix3d tetragonal_lat_matrix = cubic_lat_matrix;
        tetragonal_lat_matrix(2, 2) = 4;
        cubic_lat_ptr = std::make_unique<Lattice>(cubic_lat_matrix);
        tetragonal_lat_ptr = std::make_unique<Lattice>(tetragonal_lat_matrix);

        cu::fs::path hcp_path(cu::autotools::input_filesdir / "Mg_hcp.vasp");
        hcp_Mg_ptr = std::make_unique<Structure>(Structure::from_poscar(hcp_path));

        cu::fs::path almost_hcp_path(cu::autotools::input_filesdir / "distorted_Mg_hcp.vasp");
        almost_hcp_Mg_ptr = std::make_unique<Structure>(Structure::from_poscar(almost_hcp_path));

        rotation_90 << 0, 1, 0, -1, 0, 0, 0, 0, 1;
        translation << 0.5, 0.5, 0.5;
        eigen_vector_coordinate << 1.2, 1.3, 1.4;
        cart_op_ptr.reset(new CartOp{rotation_90, translation, false});
        coord_ptr.reset(new Coordinate{eigen_vector_coordinate});
    }

    std::unique_ptr<Lattice> cubic_lat_ptr;
    std::unique_ptr<Lattice> tetragonal_lat_ptr;
    std::unique_ptr<Structure> hcp_Mg_ptr;
    std::unique_ptr<Structure> almost_hcp_Mg_ptr;
    std::unique_ptr<CartOp> cart_op_ptr;
    std::unique_ptr<Coordinate> coord_ptr;
    Eigen::Matrix3d rotation_90;
    Eigen::Vector3d translation;
    Eigen::Vector3d eigen_vector_coordinate;
    double tol = 1e-5;
};

TEST_F(CrystalGroupTest, PointGroupSize)
{
    std::vector<cu::sym::CartOp> point_group = cu::xtal::make_point_group(primitive_fcc_Ni_ptr->lattice(), tol);
    EXPECT_EQ(point_group.size(), 48);
}

TEST_F(CrystalGroupTest, FactorGroupSize)
{
    std::vector<cu::sym::CartOp> factor_group = cu::xtal::make_factor_group(*primitive_fcc_Ni_ptr, tol);
    EXPECT_EQ(factor_group.size(), 48);
}

TEST_F(SymmetrizeTest, LatticeSymmetrize)
{
    std::vector<cu::sym::CartOp> cubic_point_group = cu::xtal::make_point_group(*cubic_lat_ptr, tol);
    cu::xtal::Lattice symmetrized_tetragonal_lattice = cu::xtal::symmetrize(*tetragonal_lat_ptr, cubic_point_group);
    std::vector<cu::sym::CartOp> symmetrized_tetragonal_point_group =
        cu::xtal::make_point_group(symmetrized_tetragonal_lattice, tol);
    EXPECT_EQ(cubic_point_group.size(), symmetrized_tetragonal_point_group.size());
}
TEST_F(SymmetrizeTest, StructureSymmetrize)
{
    std::vector<cu::sym::CartOp> hcp_factor_group = cu::xtal::make_factor_group(*hcp_Mg_ptr, tol);
    EXPECT_NE(cu::xtal::make_factor_group(*almost_hcp_Mg_ptr, tol).size(), hcp_factor_group.size());
    cu::xtal::Structure symmetrized_structure = cu::xtal::symmetrize(*almost_hcp_Mg_ptr, hcp_factor_group);
    std::vector<cu::sym::CartOp> sym_distorted_factor_group = cu::xtal::make_factor_group(symmetrized_structure, tol);
    EXPECT_EQ(sym_distorted_factor_group.size(), hcp_factor_group.size());
}

TEST_F(SymmetrizeTest, ApplySymOpEigenVector)
{
    using namespace cu::xtal;
    auto cart_op = *cart_op_ptr;
    Eigen::Vector3d transformed_eigen_vector = cart_op * eigen_vector_coordinate;
    EXPECT_TRUE(transformed_eigen_vector.isApprox(rotation_90 * eigen_vector_coordinate + translation));
}

TEST_F(SymmetrizeTest, ApplySymOpSite)
{
    cu::xtal::Site lithium_site(*coord_ptr, "Li");
    auto transformed_site = *cart_op_ptr * lithium_site;
    EXPECT_TRUE(transformed_site.cart().isApprox(rotation_90 * lithium_site.cart() + translation));
    EXPECT_EQ(transformed_site.label(), lithium_site.label());
}

TEST_F(SymmetrizeTest, ApplySymOpCoordinate)
{
    auto coord = *coord_ptr;
    auto transformed_coord = *cart_op_ptr * coord;
    EXPECT_TRUE(transformed_coord.cart().isApprox(rotation_90 * coord_ptr->cart() + translation));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
