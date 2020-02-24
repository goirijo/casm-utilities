#include <casmutils/xtal/coordinate.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>

class CoordinateTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
        coord0_ptr.reset(new rewrap::Coordinate(raw_coord));
        coord1_ptr.reset(new rewrap::Coordinate(raw_coord(0), raw_coord(1), raw_coord(2)));
    }

    // Use unique pointers because Coordinate has no default constructor
    std::unique_ptr<rewrap::Coordinate> coord0_ptr;
    std::unique_ptr<rewrap::Coordinate> coord1_ptr;
};

TEST_F(CoordinateTest, Construct) { EXPECT_EQ(coord0_ptr->cart(), coord1_ptr->cart()); }

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
