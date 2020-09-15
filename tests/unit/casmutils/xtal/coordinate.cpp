#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>

class CoordinateTest : public testing::Test
{
protected:
    using Coordinate = casmutils::xtal::Coordinate;
    using Lattice = casmutils::xtal::Lattice;
    using CoordinateEquals_f = casmutils::xtal::CoordinateEquals_f;
    // Use unique pointers because Coordinate has no default constructor
    std::unique_ptr<Coordinate> coord0_ptr;
    std::unique_ptr<Coordinate> coord1_ptr;
    std::unique_ptr<Lattice> fcc_lattice_ptr;
    std::unique_ptr<Coordinate> coord2_ptr;
    Eigen::Matrix3d lattice_matrix;
    Eigen::Vector3d frac_coords;
    double tol = 1e-5;

    void SetUp() override
    {

        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
        lattice_matrix << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
        coord0_ptr.reset(new Coordinate(raw_coord));
        coord1_ptr.reset(new Coordinate(raw_coord(0), raw_coord(1), raw_coord(2)));
        fcc_lattice_ptr.reset(new Lattice(lattice_matrix));
        frac_coords = lattice_matrix.inverse() * coord0_ptr->cart();
        coord2_ptr.reset(new Coordinate(Coordinate::from_fractional(frac_coords, *fcc_lattice_ptr)));
    }
};

TEST_F(CoordinateTest, Construct)
{
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(*coord0_ptr, *coord1_ptr, tol));
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(*coord1_ptr, *coord2_ptr, tol));
}

TEST_F(CoordinateTest, FracRetrieve) { EXPECT_TRUE(frac_coords.isApprox(coord2_ptr->frac(*fcc_lattice_ptr))); }

TEST_F(CoordinateTest, BringWithIn)
{
    // Bring within using the member function
    for (int i = -2; i <= 2; ++i)
    {
        for (int j = -2; j <= 2; ++j)
        {
            for (int l = 2; l <= 2; ++l)
            {
                Coordinate lattice_translation = Coordinate::from_fractional(i, j, l, *fcc_lattice_ptr);
                Coordinate translated_coordinate = *coord0_ptr + lattice_translation;
                translated_coordinate.bring_within(*fcc_lattice_ptr);
                EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(*coord0_ptr, translated_coordinate, tol));
            }
        }
    }
}

TEST_F(CoordinateTest, ConstBringWithIn)
{
    Coordinate lattice_translation = Coordinate::from_fractional(2, 3, 4, *fcc_lattice_ptr);
    const Coordinate translated_coordinate = *coord0_ptr + lattice_translation;
    Coordinate original_coord = translated_coordinate.bring_within(*fcc_lattice_ptr);

    EXPECT_FALSE(casmutils::is_equal<CoordinateEquals_f>(translated_coordinate, original_coord, tol));
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(*coord0_ptr, original_coord, tol));
}

TEST_F(CoordinateTest, WignerSeitzWithin)
{
    Coordinate already_within = Coordinate::from_fractional(0.25, 0.25, 0, *fcc_lattice_ptr);
    Coordinate ws_within = Coordinate::from_fractional(0.25, 0.25, 0, *fcc_lattice_ptr);

    already_within.bring_within_wigner_seitz(*fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(already_within, ws_within, tol));

    Coordinate far_right = Coordinate::from_fractional(0.75, 0.25, 0, *fcc_lattice_ptr);
    Coordinate far_far_right = Coordinate::from_fractional(1.75, 0.25, 0, *fcc_lattice_ptr);
    ws_within = Coordinate::from_fractional(-0.25, 0.25, 0, *fcc_lattice_ptr);

    far_right.bring_within_wigner_seitz(*fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(far_right, ws_within, tol));
    far_far_right.bring_within_wigner_seitz(*fcc_lattice_ptr);
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(far_far_right, ws_within, tol));
}

TEST_F(CoordinateTest, ConstWignerSeitzWithin)
{
    for (double x : {-0.1, 0.7, -2.2, 0.0})
    {
        for (double y : {-0.1, 0.7, -2.2, 0.0})
        {
            for (double z : {-0.1, 0.7, -2.2, 0.0})
            {
                const Coordinate const_coord = Coordinate::from_fractional(x, y, z, *fcc_lattice_ptr);
                Coordinate coord = const_coord;

                coord.bring_within_wigner_seitz(*fcc_lattice_ptr);
                EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(
                    const_coord.bring_within_wigner_seitz(*fcc_lattice_ptr), coord, tol));
            }
        }
    }
}

TEST_F(CoordinateTest, PlusOperator)
{
    Coordinate coord_sum = *coord0_ptr + *coord1_ptr;
    Coordinate summed_coord(coord0_ptr->cart() + coord1_ptr->cart());
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(coord_sum, summed_coord, tol));
}

TEST_F(CoordinateTest, PlusEqualOperator)
{
    Eigen::Vector3d sum_cart(0.2, 0.4, 0.6);
    Coordinate sum_coord(sum_cart);
    *coord0_ptr += *coord1_ptr;
    EXPECT_TRUE(casmutils::is_equal<CoordinateEquals_f>(*coord0_ptr, sum_coord, tol));
}

TEST_F(CoordinateTest, CoordinateEquals)
{
    CoordinateEquals_f coord0_equals(*coord0_ptr, tol);
    EXPECT_TRUE(coord0_equals(*coord1_ptr));
}

//
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
