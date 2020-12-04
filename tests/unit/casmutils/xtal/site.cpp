// these are required for test construction
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>
// tests the functions in this file
#include <casmutils/xtal/site.hpp>

class SiteTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector3d raw_coord(0.1, 0.2, 0.3);
        // store coordinate object for comparisons later
        // Makes lithium and nickel sites with the same position
        lithium_site_ptr.reset(new Site(raw_coord, "Li"));
        nickel_site_ptr.reset(new Site(raw_coord, "Ni"));
        // store lattice for coordinate transforms
        cubic_lattice_ptr.reset(new Lattice(4.0 * Eigen::Matrix3d::Identity()));
    }
    using Lattice = casmutils::xtal::Lattice;
    using Site = casmutils::xtal::Site;

    // Use unique pointers because Site has no default constructor
    std::unique_ptr<Site> lithium_site_ptr;
    std::unique_ptr<Site> nickel_site_ptr;

    std::unique_ptr<Lattice> cubic_lattice_ptr;
};

TEST_F(SiteTest, Construct)
{
    // checks construction and the const accessors
    // to the cartesian position of the site and the
    // identity of the site.
    EXPECT_EQ(lithium_site_ptr->cart(), nickel_site_ptr->cart());
    EXPECT_EQ(lithium_site_ptr->label(), "Li");
    EXPECT_EQ(nickel_site_ptr->label(), "Ni");
}

TEST_F(SiteTest, FracConversion)
{
    // checks the ability to get the fractional
    // coordinates of a site with respect to a lattice
    EXPECT_EQ(lithium_site_ptr->frac(*cubic_lattice_ptr), Eigen::Vector3d(0.025, 0.05, 0.075));
};

TEST_F(SiteTest, SiteEquals)
{
    // checks the ability to determine the equality of sites using Functor
    double tol = 1e-5;
    casmutils::xtal::SiteEquals_f is_equal_to_lithium_site(tol);
    EXPECT_TRUE(is_equal_to_lithium_site(*lithium_site_ptr, *lithium_site_ptr));
    EXPECT_FALSE(is_equal_to_lithium_site(*lithium_site_ptr, *nickel_site_ptr));

    casmutils::UnaryComparator_f<casmutils::xtal::SiteEquals_f> unary_site_comparator(*lithium_site_ptr, tol);
    EXPECT_TRUE(unary_site_comparator(*lithium_site_ptr));
    EXPECT_FALSE(unary_site_comparator(*nickel_site_ptr));
};

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
