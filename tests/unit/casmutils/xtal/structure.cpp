// These are classes that structure depends on
#include "../../../autotools.hh"
#include <casmutils/definitions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <gtest/gtest.h>
// This file tests the functions in:
#include <casmutils/xtal/structure.hpp>

class StructureTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Matrix3d cubic_lat_mat;
        cubic_lat_mat << 2, 0, 0, 0, 2, 0, 0, 0, 2;
        // Make two cubic lattices
        cubic_lat_ptr.reset(new casmutils::xtal::Lattice(cubic_lat_mat));
        big_cubic_lat_ptr.reset(new casmutils::xtal::Lattice(2 * cubic_lat_mat));
        // Make two Ni sites at different locations along z axis
        casmutils::xtal::Site site0(Eigen::Vector3d(0, 0, 1), "Ni");
        casmutils::xtal::Site site1(Eigen::Vector3d(0, 0, 2), "Ni");
        // Make bases out of each site
        basis0_ptr.reset(new std::vector<casmutils::xtal::Site>({site0}));
        basis1_ptr.reset(new std::vector<casmutils::xtal::Site>({site1}));
        // Make cubic structure
        cubic_Ni_struc_ptr.reset(new casmutils::xtal::Structure(*cubic_lat_ptr, *basis0_ptr));
        using CASM::CART;
    }

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<casmutils::xtal::Structure> cubic_Ni_struc_ptr;

    std::unique_ptr<casmutils::xtal::Lattice> cubic_lat_ptr;
    std::unique_ptr<casmutils::xtal::Lattice> big_cubic_lat_ptr;

    std::unique_ptr<std::vector<casmutils::xtal::Site>> basis0_ptr;
    std::unique_ptr<std::vector<casmutils::xtal::Site>> basis1_ptr;
    double tol = 1e-5;
};

TEST_F(StructureTest, ConstructandLatAccess)
{
    //  checks to see if Structure can be constructed
    // and has the lattice that was given associated with it.
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, cubic_Ni_struc_ptr->lattice(), tol));
}

TEST_F(StructureTest, ConstBasisAccess)
{
    //  checks the const accessor method to
    // the basis of the structure
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
        (*basis0_ptr)[0], cubic_Ni_struc_ptr->basis_sites()[0], tol));
}
TEST_F(StructureTest, ReadfromPOSCAR)
{
    //  checks the static function that reads
    // a POSCAR from a file and creates a casmutils::xtal::Structure
    // from the contents
    casmutils::fs::path pos_path(casmutils::autotools::input_filesdir / "simple_cubic_Ni.vasp");
    casmutils::xtal::Structure pos = casmutils::xtal::Structure::from_poscar(pos_path);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, pos.lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>((*basis0_ptr)[0], pos.basis_sites()[0], tol));
}
TEST_F(StructureTest, SetLatticeFrac)
{
    // FRACTIONAL CALL SHOULD CHANGE BASIS
    cubic_Ni_struc_ptr->set_lattice(*big_cubic_lat_ptr, casmutils::xtal::FRAC);
    // Basis should now be different (site0 transforms to site1 from test frame)
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
        (*basis1_ptr)[0], cubic_Ni_struc_ptr->basis_sites()[0], tol));
    // Lattice should be different
    EXPECT_FALSE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, cubic_Ni_struc_ptr->lattice(), tol));
}
TEST_F(StructureTest, ConstSetLatticeCart)
{
    //  checks the set lattice function
    // in one instance keeping the fractional coordinates
    // and in one instance keeping the cartesian coordinates

    // CARTESIAN CALL SHOULD NOT CHANGE BASIS
    const casmutils::xtal::Structure& cubic_Ni_strucref = *cubic_Ni_struc_ptr;
    casmutils::xtal::Structure new_struc = cubic_Ni_strucref.set_lattice(*big_cubic_lat_ptr, casmutils::xtal::CART);
    // Const call to set_lattice doesn't mutate original
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::SiteEquals_f>((*basis0_ptr)[0], cubic_Ni_strucref.basis_sites()[0], tol));
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, cubic_Ni_strucref.lattice(), tol));
    // Const call to set_lattice does return new structure
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>((*basis0_ptr)[0], new_struc.basis_sites()[0], tol));
    EXPECT_FALSE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, new_struc.lattice(), tol));
}
TEST_F(StructureTest, Within)
{
    // checks to see if all basis sites can be moved within the bounding box of the lattice
    casmutils::xtal::Site outside_site(Eigen::Vector3d(0, 0, -1), "Ni");
    casmutils::xtal::Structure outside_structure(*cubic_lat_ptr, {outside_site});
    outside_structure.within();
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*cubic_lat_ptr, outside_structure.lattice(), tol));
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::SiteEquals_f>((*basis0_ptr)[0], outside_structure.basis_sites()[0], tol));
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
