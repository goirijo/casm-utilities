#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>

class CoordinateTest : public testing::Test
{
protected:
    using Lattice = casmutils::xtal::Lattice;
    using CoordinateEquals_f = casmutils::xtal::CoordinateEquals_f;
    // Use unique pointers because Coordinate has no default constructor
    Eigen::Vector3d coord0;
    std::unique_ptr<Lattice> fcc_lattice_ptr;
    Eigen::Matrix3d lattice_matrix;
    Eigen::Vector3d frac_coords;
    double tol = 1e-5;

    void SetUp() override
    {
        coord0 << 0.1, 0.2, 0.3;
        lattice_matrix << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
        fcc_lattice_ptr.reset(new Lattice(lattice_matrix));
        frac_coords = lattice_matrix.inverse() * coord0;
    }
};

TEST_F(CoordinateTest, CartRetrieve)
{
    Eigen::Vector3d new_cart_coords = casmutils::xtal::fractional_to_cartesian(frac_coords, *fcc_lattice_ptr);
    EXPECT_TRUE(coord0.isApprox(new_cart_coords, tol));
}

TEST_F(CoordinateTest, FracRetrieve)
{
    Eigen::Vector3d new_frac_coords = casmutils::xtal::cartesian_to_fractional(coord0, *fcc_lattice_ptr);
    EXPECT_TRUE(frac_coords.isApprox(new_frac_coords, tol));
}

TEST_F(CoordinateTest, BringWithIn)
{
    // Bring within using the member function
    for (int i = -2; i <= 2; ++i)
    {
        for (int j = -2; j <= 2; ++j)
        {
            for (int l = 2; l <= 2; ++l)
            {
                Eigen::Vector3d lattice_translation =
                    casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(i, j, l), *fcc_lattice_ptr);
                Eigen::Vector3d translated_coordinate = coord0 + lattice_translation;
                Eigen::Vector3d withined_coords =
                    casmutils::xtal::bring_within_lattice(translated_coordinate, *fcc_lattice_ptr);
                EXPECT_TRUE(coord0.isApprox(withined_coords));
            }
        }
    }
}

TEST_F(CoordinateTest, WignerSeitzWithin)
{
    auto already_within = casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(0.25, 0.25, 0), *fcc_lattice_ptr);
    auto ws_within = casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(0.25, 0.25, 0), *fcc_lattice_ptr);

    auto new_within = casmutils::xtal::bring_within_wigner_seitz(already_within, *fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(new_within, ws_within, tol));

    auto far_right = casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(0.75, 0.25, 0), *fcc_lattice_ptr);
    auto far_far_right = casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(1.75, 0.25, 0), *fcc_lattice_ptr);
    ws_within = casmutils::xtal::fractional_to_cartesian(Eigen::Vector3d(-0.25, 0.25, 0), *fcc_lattice_ptr);

    auto new_far_right = casmutils::xtal::bring_within_wigner_seitz(far_right, *fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(new_far_right, ws_within, tol));

    auto new_far_far_right = casmutils::xtal::bring_within_wigner_seitz(far_far_right, *fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(new_far_far_right, ws_within, tol));
}

TEST_F(CoordinateTest, CoordinateEquals)
{
    CoordinateEquals_f coord0_equals(tol);
    EXPECT_TRUE(coord0_equals(coord0, coord0));

    casmutils::UnaryComparator_f<CoordinateEquals_f> unary_coord_compare(coord0, tol);
    EXPECT_TRUE(unary_coord_compare(coord0));
}

//
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
