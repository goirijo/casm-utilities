#include "../../../autotools.hh"
#include "casmutils/sym/cartesian.hpp"
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
