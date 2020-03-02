// These are classes that structure depends on
#include "../../autotools.hh"
#include <casmutils/definitions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <gtest/gtest.h>
// This file tests the functions in:
#include <casmutils/xtal/structure_tools.hpp>

class StructureToolsTest : public testing::Test
{
protected:
    using Structure = casmutils::xtal::Structure;
    void SetUp() override
    {
        // Paths to testing poscars
        casmutils::fs::path cubic_path(casmutils::autotools::input_filesdir / "simple_cubic_Ni.vasp");
        casmutils::fs::path conventional_path(casmutils::autotools::input_filesdir / "conventional_fcc_Ni.vasp");
        casmutils::fs::path nonniggli_conventional_path(casmutils::autotools::input_filesdir /
                                                        "nonniggli_conventional_fcc_Ni.vasp");
        casmutils::fs::path primitive_path(casmutils::autotools::input_filesdir / "primitive_fcc_Ni.vasp");
        casmutils::fs::path deformed_conventional_path(casmutils::autotools::input_filesdir /
                                                       "deformed_conventional_fcc_Ni.vasp");
        casmutils::fs::path strained_conventional_path(casmutils::autotools::input_filesdir /
                                                       "strained_conventional_fcc_Ni.vasp");
        // test deformation
        deformation << 1.0, 0.02, 0.02, 0.0, 1.0, 0.0, 0.0, 0.0, 1.05;
        // test GL strain
        Eigen::Matrix3d GL_tensor = (deformation * deformation - Eigen::Matrix3d::Identity()) / 2;
        GL_unrolled_strain.resize(6);
        GL_unrolled_strain << GL_tensor(0, 0), GL_tensor(1, 1), GL_tensor(2, 2), GL_tensor(1, 2), GL_tensor(0, 2),
            GL_tensor(0, 1);
        // load the poscars
        cubic_Ni_struc_ptr = std::make_unique<Structure>(Structure::from_poscar(cubic_path));
        conventional_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(conventional_path));
        nonniggli_conventional_fcc_Ni_ptr =
            std::make_unique<Structure>(Structure::from_poscar(nonniggli_conventional_path));
        primitive_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(primitive_path));
        deformed_conventional_fcc_Ni_ptr =
            std::make_unique<Structure>(Structure::from_poscar(deformed_conventional_path));
        strained_conventional_fcc_Ni_ptr =
            std::make_unique<Structure>(Structure::from_poscar(strained_conventional_path));
    }

    // returns true if ref basis is exactly the same as test basis (permutation not allowed) false otherwise
    bool cartesian_basis_is_equal(const std::vector<casmutils::xtal::Site>& ref_basis,
                                  const std::vector<casmutils::xtal::Site>& test_basis)
    {
        if (ref_basis.size() != test_basis.size())
        {
            return false;
        }
        for (int i = 0; i < ref_basis.size(); i++)
        {
            if (!casmutils::is_equal<casmutils::xtal::SiteEquals_f>(ref_basis[i], test_basis[i], tol))
            {
                return false;
            }
        }
        return true;
    }

    // returns true if fractional ref basis is exactly the same as fractional test basis (permutation not allowed) false
    // otherwise
    bool fractional_basis_is_equal(const casmutils::xtal::Lattice& ref_lattice,
                                   const std::vector<casmutils::xtal::Site>& ref_basis,
                                   const casmutils::xtal::Lattice& test_lattice,
                                   const std::vector<casmutils::xtal::Site>& test_basis)
    {
        if (ref_basis.size() != test_basis.size())
        {
            return false;
        }
        for (int i = 0; i < ref_basis.size(); i++)
        {
            if (!ref_basis[i].frac(ref_lattice).isApprox(test_basis[i].frac(test_lattice), tol))
            {
                return false;
            }
        }
        return true;
    }
    // returns true if ref basis is the same as test basis (permutation allowed) false otherwise
    bool cartesian_basis_is_equal_with_permutation(const std::vector<casmutils::xtal::Site>& ref_basis,
                                                   const std::vector<casmutils::xtal::Site>& test_basis)
    {
        if (ref_basis.size() != test_basis.size())
        {
            return false;
        }
        for (int i = 0; i < ref_basis.size(); i++)
        {
            // construct an equals predicate
            casmutils::xtal::SiteEquals_f is_equal_to_site_i(ref_basis[i], tol);
            // search in constructed basis for site equivalent to conventional fcc site i
            if (std::find_if(test_basis.begin(), test_basis.end(), is_equal_to_site_i) == test_basis.end())
            {
                // if hit end iterator return false
                return false;
            }
        }
        return true;
    }

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<Structure> cubic_Ni_struc_ptr;
    std::unique_ptr<Structure> conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> nonniggli_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> deformed_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> strained_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> primitive_fcc_Ni_ptr;

    double tol = 1e-5;

    Eigen::Matrix3d deformation;
    Eigen::VectorXd GL_unrolled_strain;
};

TEST_F(StructureToolsTest, WritePOSCAR)
{
    // checks to see if writing a POSCAR from a structure
    // can be read as the same structure
    casmutils::fs::path write_path(casmutils::autotools::input_filesdir / "simple_cubic_Ni_copy.vasp");
    casmutils::xtal::write_poscar(*cubic_Ni_struc_ptr, write_path);
    const Structure pos = Structure::from_poscar(write_path);
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(cubic_Ni_struc_ptr->lattice(), pos.lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(cubic_Ni_struc_ptr->basis_sites(), pos.basis_sites()));
}
TEST_F(StructureToolsTest, MakePrimitive)
{
    // checks to see if conventional fcc gets reduced to a primitive fcc
    const Structure constructed_primitive = casmutils::xtal::make_primitive(*conventional_fcc_Ni_ptr);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(primitive_fcc_Ni_ptr->lattice(),
                                                                      constructed_primitive.lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(primitive_fcc_Ni_ptr->basis_sites(), constructed_primitive.basis_sites()));
}
TEST_F(StructureToolsTest, MakeSuperstructure)
{
    // checks to see if primitive fcc transforms into conventional cell correctly
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    const Structure constructed_superstructure =
        casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, transf_mat);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(conventional_fcc_Ni_ptr->lattice(),
                                                                      constructed_superstructure.lattice(), tol));
    const auto& constructed_basis = constructed_superstructure.basis_sites();
    // need to compare basis differently here because it may be permuted
    EXPECT_TRUE(cartesian_basis_is_equal_with_permutation(conventional_fcc_Ni_ptr->basis_sites(),
                                                          constructed_superstructure.basis_sites()));
}
TEST_F(StructureToolsTest, MakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    casmutils::xtal::make_niggli(nonniggli_conventional_fcc_Ni_ptr.get());
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(
        conventional_fcc_Ni_ptr->lattice(), nonniggli_conventional_fcc_Ni_ptr->lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(conventional_fcc_Ni_ptr->basis_sites(),
                                         nonniggli_conventional_fcc_Ni_ptr->basis_sites()));
}
TEST_F(StructureToolsTest, ConstMakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    const Structure niggli = casmutils::xtal::make_niggli(*nonniggli_conventional_fcc_Ni_ptr);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(conventional_fcc_Ni_ptr->lattice(),
                                                                      niggli.lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(conventional_fcc_Ni_ptr->basis_sites(), niggli.basis_sites()));
}
TEST_F(StructureToolsTest, ApplyDeformation)
{
    // checks the application of a deformation matrix on a structure
    casmutils::xtal::apply_deformation(conventional_fcc_Ni_ptr.get(), deformation);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(deformed_conventional_fcc_Ni_ptr->lattice(),
                                                                      conventional_fcc_Ni_ptr->lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(deformed_conventional_fcc_Ni_ptr->basis_sites(),
                                         conventional_fcc_Ni_ptr->basis_sites()));
}
TEST_F(StructureToolsTest, ConstApplyDeformation)
{
    // checks the application of a deformation matrix on a structure
    const Structure deformed_fcc_Ni = casmutils::xtal::apply_deformation(*conventional_fcc_Ni_ptr, deformation);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(deformed_conventional_fcc_Ni_ptr->lattice(),
                                                                      deformed_fcc_Ni.lattice(), tol));
    EXPECT_TRUE(
        cartesian_basis_is_equal(deformed_conventional_fcc_Ni_ptr->basis_sites(), deformed_fcc_Ni.basis_sites()));
    // cartesian basis should change on deformation
    EXPECT_FALSE(cartesian_basis_is_equal(conventional_fcc_Ni_ptr->basis_sites(), deformed_fcc_Ni.basis_sites()));
    EXPECT_TRUE(fractional_basis_is_equal(conventional_fcc_Ni_ptr->lattice(), conventional_fcc_Ni_ptr->basis_sites(),
                                          deformed_fcc_Ni.lattice(), deformed_fcc_Ni.basis_sites()));
}
TEST_F(StructureToolsTest, ApplyStrain)
{
    // checks the application of a GL strain on a structure
    casmutils::xtal::apply_strain(conventional_fcc_Ni_ptr.get(), GL_unrolled_strain, "GL");
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(strained_conventional_fcc_Ni_ptr->lattice(),
                                                                      conventional_fcc_Ni_ptr->lattice(), tol));
    EXPECT_TRUE(cartesian_basis_is_equal(strained_conventional_fcc_Ni_ptr->basis_sites(),
                                         conventional_fcc_Ni_ptr->basis_sites()));
}
TEST_F(StructureToolsTest, ConstApplyStrain)
{
    // checks the application of a GL strain on a structure
    const Structure strained_fcc_Ni = casmutils::xtal::apply_strain(*conventional_fcc_Ni_ptr, GL_unrolled_strain, "GL");
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(strained_conventional_fcc_Ni_ptr->lattice(),
                                                                      strained_fcc_Ni.lattice(), tol));
    EXPECT_TRUE(
        cartesian_basis_is_equal(strained_conventional_fcc_Ni_ptr->basis_sites(), strained_fcc_Ni.basis_sites()));
    // cartesian basis should change on straining
    EXPECT_FALSE(cartesian_basis_is_equal(conventional_fcc_Ni_ptr->basis_sites(), strained_fcc_Ni.basis_sites()));
    // fractional basis should remain unchanged
    EXPECT_TRUE(fractional_basis_is_equal(conventional_fcc_Ni_ptr->lattice(), conventional_fcc_Ni_ptr->basis_sites(),
                                          strained_fcc_Ni.lattice(), strained_fcc_Ni.basis_sites()));
}
TEST_F(StructureToolsTest, MakeSuperstructuresofVol)
{
    // create transformation matrices for volume 2 and 3
    Eigen::Matrix3i vol2_var0, vol2_var1, vol3_var0, vol3_var1, vol3_var2;
    vol2_var0 << 1, 0, 1, 0, 1, 1, -1, -1, 0;
    vol2_var1 << 0, -1, -1, 1, 0, 1, 0, 1, -1;
    vol3_var0 << -1, 1, 1, 0, -1, 1, 1, 0, 1;
    vol3_var1 << -1, 1, 0, 0, -1, 2, 1, 1, -1;
    vol3_var2 << -1, 0, 2, 0, 1, -2, 1, 0, 1;
    // create structures
    Structure struc_vol2_0 =
        casmutils::xtal::make_niggli(casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, vol2_var0));
    Structure struc_vol2_1 =
        casmutils::xtal::make_niggli(casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, vol2_var1));
    Structure struc_vol3_0 =
        casmutils::xtal::make_niggli(casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, vol3_var0));
    Structure struc_vol3_1 =
        casmutils::xtal::make_niggli(casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, vol3_var1));
    Structure struc_vol3_2 =
        casmutils::xtal::make_niggli(casmutils::xtal::make_super_structure(*primitive_fcc_Ni_ptr, vol3_var2));
    // create vol2 and vol3 all at same time
    std::vector<Structure> vol2_superstrucs = casmutils::xtal::make_superstructures_of_volume(*primitive_fcc_Ni_ptr, 2);
    std::vector<Structure> vol3_superstrucs = casmutils::xtal::make_superstructures_of_volume(*primitive_fcc_Ni_ptr, 3);
    // These equalities are sensitive to similarity transforms as well as the ordering the in vector
    // More robust equality checking methods could ensure that this test lasts
    // check lattices
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(struc_vol2_0.lattice(),
                                                                      vol2_superstrucs[0].lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(struc_vol2_1.lattice(),
                                                                      vol2_superstrucs[1].lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(struc_vol3_0.lattice(),
                                                                      vol3_superstrucs[0].lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(struc_vol3_1.lattice(),
                                                                      vol3_superstrucs[1].lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(struc_vol3_2.lattice(),
                                                                      vol3_superstrucs[2].lattice(), tol));
    // check bases
    EXPECT_TRUE(
        cartesian_basis_is_equal_with_permutation(struc_vol2_0.basis_sites(), vol2_superstrucs[0].basis_sites()));
    EXPECT_TRUE(
        cartesian_basis_is_equal_with_permutation(struc_vol2_1.basis_sites(), vol2_superstrucs[1].basis_sites()));
    EXPECT_TRUE(
        cartesian_basis_is_equal_with_permutation(struc_vol3_0.basis_sites(), vol3_superstrucs[0].basis_sites()));
    EXPECT_TRUE(
        cartesian_basis_is_equal_with_permutation(struc_vol3_1.basis_sites(), vol3_superstrucs[1].basis_sites()));
    EXPECT_TRUE(
        cartesian_basis_is_equal_with_permutation(struc_vol3_2.basis_sites(), vol3_superstrucs[2].basis_sites()));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
