#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <utility>

class CoordinateTest : public testing::Test{
protected:
    void SetUp() override{
        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
        lattice_matrix << 0,0.5,0.5,0.5,0,0.5,0.5,0.5,0;
        coord0_ptr.reset(new casmutils::xtal::Coordinate(raw_coord));
        coord1_ptr.reset(new casmutils::xtal::Coordinate(raw_coord(0), raw_coord(1), raw_coord(2)));
        lattice_ptr.reset(new casmutils::xtal::Lattice(lattice_matrix));
        frac_coords = lattice_matrix.inverse()*coord0_ptr->cart(); 
    }
    // Use unique pointers because Coordinate has no default constructor
    std::unique_ptr<casmutils::xtal::Coordinate> coord0_ptr;
    std::unique_ptr<casmutils::xtal::Coordinate> coord1_ptr;
    std::unique_ptr<casmutils::xtal::Lattice> lattice_ptr;
    Eigen::Matrix3d lattice_matrix;
    Eigen::Vector3d frac_coords;
};

TEST_F(CoordinateTest, Construct) { EXPECT_EQ(coord0_ptr->cart(), coord1_ptr->cart()); }

TEST_F(CoordinateTest, FracRetrieve){ 
    casmutils::xtal::Coordinate coord2 = casmutils::xtal::Coordinate::from_fractional(frac_coords,*lattice_ptr);
    EXPECT_TRUE(frac_coords.isApprox(coord2.frac(*lattice_ptr)));
}

TEST_F(CoordinateTest, BringWithIn){
    //Bring within using the member function
    casmutils::xtal::Coordinate coord2 = casmutils::xtal::Coordinate::from_fractional(frac_coords,*lattice_ptr);
    coord2.bring_within(*lattice_ptr);  
    //Bring within using formuation
    for (int i=0; i<3; ++i){
        if (frac_coords(i)>=1 || frac_coords(i)<0){
            frac_coords(i) = frac_coords(i)-floor(frac_coords(i)); 
        }
    }
    
    EXPECT_TRUE(coord2.frac(*lattice_ptr).isApprox(frac_coords));
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
