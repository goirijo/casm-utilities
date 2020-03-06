#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>

class CoordinateTest : public testing::Test{
protected:
    void SetUp() override{
        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
        Eigen::Matrix3d lattice;
        lattice << 0,0.5,0.5,
                0.5,0,0.5,
                0.5,0.5,0;
        
        coord0_ptr.reset(new rewrap::Coordinate(raw_coord));
        coord1_ptr.reset(new rewrap::Coordinate(raw_coord(0), raw_coord(1), raw_coord(2)));
        lattice_ptr.reset(new rewrap::Lattice(lattice));
    }
    // Use unique pointers because Coordinate has no default constructor
    std::unique_ptr<rewrap::Coordinate> coord0_ptr;
    std::unique_ptr<rewrap::Coordinate> coord1_ptr;
    std::unique_ptr<rewrap::Lattice> lattice_ptr;
};

TEST_F(CoordinateTest, Construct) { EXPECT_EQ(coord0_ptr->cart(), coord1_ptr->cart()); }

TEST_F(CoordinateTest, FracRetrieve){
    rewrap::Lattice &lattice_rewrap = *lattice_ptr;
    Eigen::Matrix3d lattice;
    lattice << 0,0.5,0.5,0.5,0,0.5,0.5,0.5,0;
    Eigen::Vector3d frac_coords = lattice.inverse()*coord0_ptr->cart();
  
    rewrap::Coordinate coords_rewrap = rewrap::Coordinate::from_fractional(frac_coords,lattice_rewrap);     
    EXPECT_TRUE(frac_coords.isApprox(coords_rewrap.frac(lattice_rewrap)));
}

TEST_F(CoordinateTest, BringWithIn){
    Eigen::Matrix3d lattice;
    lattice << 0,0.5,0.5,0.5,0,0.5,0.5,0.5,0;
    Eigen::Vector3d frac_coords = lattice.inverse()*coord0_ptr->cart();
    //Bring within using the member function
    rewrap::Coordinate coords_rewrap = rewrap::Coordinate::from_fractional(frac_coords,*lattice_ptr);
    coords_rewrap.bring_within(*lattice_ptr);  
    //Bring within using formuation
    for (int i=0; i<3; ++i){
        if (frac_coords(i)>=1 || frac_coords(i)<0){
            frac_coords(i) = frac_coords(i)-floor(frac_coords(i)); 
        }
    }
    
    EXPECT_TRUE(coords_rewrap.frac(*lattice_ptr).isApprox(frac_coords));
}
TEST_F(CoordinateTest, PlusEqualOperator){
    Eigen::Vector3d coord_to_compare(0.2,0.4,0.6);
    *coord0_ptr += *coord1_ptr;
    EXPECT_EQ(coord0_ptr->cart(),coord_to_compare);
}

TEST_F(CoordinateTest, EqualOperator){
    EXPECT_TRUE(*coord0_ptr==*coord1_ptr);
}
//
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
