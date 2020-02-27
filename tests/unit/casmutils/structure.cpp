// These are classes that structure depends on
#include <casmutils/definitions.hpp>
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
        cubic_lat_ptr.reset(new casmutils::xtal::Lattice(cubic_lat_mat));
        double_cubic_lat_ptr.reset(new casmutils::xtal::Lattice(2 * cubic_lat_mat));
        lattice_equal_ptr.reset(new casmutils::xtal::LatticeEquals_f(*cubic_lat_ptr, tol));
        casmutils::xtal::Site site0(casmutils::xtal::Coordinate(Eigen::Vector3d(0, 0, 1)), "Ni");
        casmutils::xtal::Site site1(casmutils::xtal::Coordinate(Eigen::Vector3d(0, 0, 2)), "Ni");
        basis0_ptr.reset(new std::vector<casmutils::xtal::Site>({site0}));
        basis1_ptr.reset(new std::vector<casmutils::xtal::Site>({site1}));
        struc0_ptr.reset(new casmutils::xtal::Structure(*cubic_lat_ptr, *basis0_ptr));
        site0_equal_ptr.reset(new casmutils::xtal::SiteEquals_f((*basis0_ptr)[0], tol));
        site1_equal_ptr.reset(new casmutils::xtal::SiteEquals_f((*basis1_ptr)[0], tol));
        // struc1_ptr.reset(new casmutils::xtal::Structure(*cubic_lat_ptr,*basis0_ptr));
        using CASM::CART;
    }

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<casmutils::xtal::Structure> struc0_ptr;
    // std::unique_ptr<casmutils::xtal::Structure> struc1_ptr;
    std::unique_ptr<casmutils::xtal::Lattice> cubic_lat_ptr;
    std::unique_ptr<casmutils::xtal::Lattice> double_cubic_lat_ptr;
    std::unique_ptr<casmutils::xtal::LatticeEquals_f> lattice_equal_ptr;
    std::unique_ptr<std::vector<casmutils::xtal::Site>> basis0_ptr;
    std::unique_ptr<std::vector<casmutils::xtal::Site>> basis1_ptr;
    std::unique_ptr<casmutils::xtal::SiteEquals_f> site0_equal_ptr;
    std::unique_ptr<casmutils::xtal::SiteEquals_f> site1_equal_ptr;
    double tol = 1e-5;
};

TEST_F(StructureTest, ConstructandLatAccess)
{
    // This test checks to see if Structure can be constructed
    // and has the lattice that was given associated with it.
    EXPECT_TRUE((*lattice_equal_ptr)(struc0_ptr->lattice()));
}

TEST_F(StructureTest, ConstBasisAccess)
{
    // This test checks the const accessor method to
    // the basis of the structure
    EXPECT_TRUE((*site0_equal_ptr)(struc0_ptr->basis_sites()[0]));
}
TEST_F(StructureTest,ReadfromPOSCAR){
	// This test checks the static function that reads 
	// a POSCAR from a file and creates a casmutils::xtal::Structure
	// from the contents
	casmutils::fs::path pos_path("input_files/simple_cubic_Ni.vasp");
	casmutils::xtal::Structure pos=casmutils::xtal::Structure::from_poscar(pos_path.string());
    EXPECT_TRUE((*lattice_equal_ptr)(pos.lattice()));
    EXPECT_TRUE((*site0_equal_ptr)(pos.basis_sites()[0]));
}
TEST_F(StructureTest, SetLattice)
{
    // This test checks the set lattice function
    // in one instance keeping the fractional coordinates
    // and in one instance keeping the cartesian coordinates

    // CARTESIAN CALL SHOULD NOT CHANGE BASIS
    const casmutils::xtal::Structure& struc0ref = *struc0_ptr;
    casmutils::xtal::Structure new_struc = struc0ref.set_lattice(*double_cubic_lat_ptr, casmutils::xtal::CART);
    // Const call to set_lattice doesn't mutate original
    EXPECT_TRUE((*site0_equal_ptr)(struc0ref.basis_sites()[0]));
    EXPECT_TRUE((*lattice_equal_ptr)(struc0ref.lattice()));
    // Const call to set_lattice does return new structure
    EXPECT_TRUE((*site0_equal_ptr)(new_struc.basis_sites()[0]));
    EXPECT_TRUE(!(*lattice_equal_ptr)(new_struc.lattice()));
    /// FRACTIONAL CALL SHOULD CHANGE BASIS
    struc0_ptr->set_lattice(*double_cubic_lat_ptr, casmutils::xtal::FRAC);
    /// Basis should now be different (site0 transforms to site1 from test frame)
    EXPECT_TRUE((*site1_equal_ptr)(struc0_ptr->basis_sites()[0]));
    /// Lattice should be different
    EXPECT_TRUE(!(*lattice_equal_ptr)(struc0_ptr->lattice()));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
