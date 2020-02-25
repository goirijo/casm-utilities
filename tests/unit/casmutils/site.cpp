//these are required for test construction
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
//tests the functions in this file
#include <casmutils/xtal/site.hpp>
class SiteTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
		rewrap::Coordinate init_position(raw_coord);
		//store coordinate object for comparisons later
		init_pos_ptr.reset(new rewrap::Coordinate(raw_coord));

        lithium_site_ptr.reset(new rewrap::Site(init_position,"Li"));
        nickel_site_ptr.reset(new rewrap::Site(init_position,"Ni"));
		//store lattice for coordinate transforms
		cubic_lattice_ptr.reset(new rewrap::Lattice(4.0*Eigen::Matrix3d::Identity()));

    }

    // Use unique pointers because Site has no default constructor
    std::unique_ptr<rewrap::Site> lithium_site_ptr;
    std::unique_ptr<rewrap::Site> nickel_site_ptr;

    std::unique_ptr<rewrap::Coordinate> init_pos_ptr;
    std::unique_ptr<rewrap::Lattice> cubic_lattice_ptr;
};


TEST_F(SiteTest, Construct) { 
	//This test checks construction and the const accessors
	// to the cartesian position of the site and the 
	// identity of the site.	
	EXPECT_EQ(lithium_site_ptr->cart(), nickel_site_ptr->cart()); 
	EXPECT_EQ(lithium_site_ptr->label(),"Li");
	EXPECT_EQ(nickel_site_ptr->label(),"Ni");
}

TEST_F(SiteTest, CoordCast) {
	//This test checks the ability to cast a Site to a Coordinate
	EXPECT_EQ(rewrap::Coordinate(*lithium_site_ptr).cart(),init_pos_ptr->cart());
};

TEST_F(SiteTest, FracConversion) {
	//This test checks the ability to get the fractional
	// coordinates of a site with respect to a lattice
	EXPECT_EQ(lithium_site_ptr->frac(*cubic_lattice_ptr),Eigen::Vector3d(0.025,0.05,0.075));
};

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
