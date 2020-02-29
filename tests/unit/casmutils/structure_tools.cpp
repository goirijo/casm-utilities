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
        GL_unrolled_strain = {GL_tensor(0, 0), GL_tensor(1, 1), GL_tensor(2, 2),
                              GL_tensor(1, 2), GL_tensor(0, 2), GL_tensor(0, 1)};

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

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<Structure> cubic_Ni_struc_ptr;
    std::unique_ptr<Structure> conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> nonniggli_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> deformed_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> strained_conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> primitive_fcc_Ni_ptr;

    double tol = 1e-5;

    Eigen::Matrix3d deformation;
    std::vector<double> GL_unrolled_strain;
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
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(cubic_Ni_struc_ptr->basis_sites()[0],
                                                                   pos.basis_sites()[0], tol));
}
TEST_F(StructureToolsTest, MakePrimitive)
{
    // checks to see if conventional fcc gets reduced to a primitive fcc
    const Structure constructed_primitive = casmutils::xtal::make_primitive(*conventional_fcc_Ni_ptr);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(primitive_fcc_Ni_ptr->lattice(),
                                                                      constructed_primitive.lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(primitive_fcc_Ni_ptr->basis_sites()[0],
                                                                   constructed_primitive.basis_sites()[0], tol));
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
    for (int i = 0; i < 4; i++)
    {
        // construct an equals predicate
        casmutils::xtal::SiteEquals_f is_equal_to_site_i(conventional_fcc_Ni_ptr->basis_sites()[i], tol);
        // search in constructed basis for site equivalent to conventional fcc site i
        EXPECT_NE(std::find_if(constructed_basis.begin(), constructed_basis.end(), is_equal_to_site_i),
                  constructed_basis.end());
    }
}
TEST_F(StructureToolsTest, MakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    casmutils::xtal::make_niggli(nonniggli_conventional_fcc_Ni_ptr.get());
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(
        conventional_fcc_Ni_ptr->lattice(), nonniggli_conventional_fcc_Ni_ptr->lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
            conventional_fcc_Ni_ptr->basis_sites()[i], nonniggli_conventional_fcc_Ni_ptr->basis_sites()[i], tol));
    }
}
TEST_F(StructureToolsTest, ConstMakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    const Structure niggli = casmutils::xtal::make_niggli(*nonniggli_conventional_fcc_Ni_ptr);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(conventional_fcc_Ni_ptr->lattice(),
                                                                      niggli.lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(conventional_fcc_Ni_ptr->basis_sites()[i],
                                                                       niggli.basis_sites()[i], tol));
    }
}
TEST_F(StructureToolsTest, ApplyDeformation)
{
    // checks the application of a deformation matrix on a structure
    casmutils::xtal::apply_deformation(conventional_fcc_Ni_ptr.get(), deformation);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(deformed_conventional_fcc_Ni_ptr->lattice(),
                                                                      conventional_fcc_Ni_ptr->lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
            deformed_conventional_fcc_Ni_ptr->basis_sites()[i], conventional_fcc_Ni_ptr->basis_sites()[i], tol));
    }
}
TEST_F(StructureToolsTest, ConstApplyDeformation)
{
    // checks the application of a deformation matrix on a structure
    const Structure deformed_fcc_Ni = casmutils::xtal::apply_deformation(*conventional_fcc_Ni_ptr, deformation);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(deformed_conventional_fcc_Ni_ptr->lattice(),
                                                                      deformed_fcc_Ni.lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
            deformed_conventional_fcc_Ni_ptr->basis_sites()[i], deformed_fcc_Ni.basis_sites()[i], tol));
        // cartesian location of sites should have changed
        if (i != 0)
        {
            // every point but the origin changes its cartesian position
            EXPECT_FALSE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(conventional_fcc_Ni_ptr->basis_sites()[i],
                                                                            deformed_fcc_Ni.basis_sites()[i], tol));
        }
        else
        {
            // origin point is index 0
            EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(conventional_fcc_Ni_ptr->basis_sites()[i],
                                                                           deformed_fcc_Ni.basis_sites()[i], tol));
        }
        // fractional sites should not change
        bool isequal = conventional_fcc_Ni_ptr->basis_sites()[i]
                           .frac(conventional_fcc_Ni_ptr->lattice())
                           .isApprox(deformed_fcc_Ni.basis_sites()[i].frac(deformed_fcc_Ni.lattice()), tol);
        EXPECT_TRUE(isequal);
    }
}
TEST_F(StructureToolsTest, ApplyStrain)
{
    // checks the application of a GL strain on a structure
    Eigen::VectorXd casted =
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(GL_unrolled_strain.data(), GL_unrolled_strain.size());
    casmutils::xtal::apply_strain(conventional_fcc_Ni_ptr.get(), casted, "GL");
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(strained_conventional_fcc_Ni_ptr->lattice(),
                                                                      conventional_fcc_Ni_ptr->lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
            strained_conventional_fcc_Ni_ptr->basis_sites()[i], conventional_fcc_Ni_ptr->basis_sites()[i], tol));
    }
}
TEST_F(StructureToolsTest, ConstApplyStrain)
{
    // checks the application of a GL strain on a structure
    Eigen::VectorXd casted =
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(GL_unrolled_strain.data(), GL_unrolled_strain.size());
    const Structure strained_fcc_Ni = casmutils::xtal::apply_strain(*conventional_fcc_Ni_ptr, casted, "GL");
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(strained_conventional_fcc_Ni_ptr->lattice(),
                                                                      strained_fcc_Ni.lattice(), tol));
    for (int i = 0; i < 4; i++)
    {
        EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(
            strained_conventional_fcc_Ni_ptr->basis_sites()[i], strained_fcc_Ni.basis_sites()[i], tol));
        // cartesian location of sites should have changed
        if (i != 0)
        {
            // every point but the origin changes its cartesian position
            EXPECT_FALSE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(conventional_fcc_Ni_ptr->basis_sites()[i],
                                                                            strained_fcc_Ni.basis_sites()[i], tol));
        }
        else
        {
            // origin point is index 0
            EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(conventional_fcc_Ni_ptr->basis_sites()[i],
                                                                           strained_fcc_Ni.basis_sites()[i], tol));
        }
        // fractional sites should not change
        bool isequal = conventional_fcc_Ni_ptr->basis_sites()[i]
                           .frac(conventional_fcc_Ni_ptr->lattice())
                           .isApprox(strained_fcc_Ni.basis_sites()[i].frac(strained_fcc_Ni.lattice()), tol);
        EXPECT_TRUE(isequal);
    }
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
