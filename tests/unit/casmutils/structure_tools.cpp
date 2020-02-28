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
    void SetUp() override
    {

        casmutils::fs::path pos_path(casmutils::autotools::input_filesdir / "simple_cubic_Ni.vasp");
        cubic_Ni_struc_ptr.reset(new casmutils::xtal::Structure(casmutils::xtal::Structure::from_poscar(pos_path)));
    }

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<casmutils::xtal::Structure> cubic_Ni_struc_ptr;
    double tol = 1e-5;
};

TEST_F(StructureToolsTest, WritePOSCAR)
{
    // checks to see if writing a POSCAR from a structure
    // can be read as the same structure
    casmutils::fs::path write_path(casmutils::autotools::input_filesdir / "simple_cubic_Ni_copy.vasp");
    casmutils::xtal::write_poscar(*cubic_Ni_struc_ptr, write_path);
    const casmutils::xtal::Structure pos = casmutils::xtal::Structure::from_poscar(write_path);
    EXPECT_TRUE(
        casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(cubic_Ni_struc_ptr->lattice(), pos.lattice(), tol));
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::SiteEquals_f>(cubic_Ni_struc_ptr->basis_sites()[0],
                                                                   pos.basis_sites()[0], tol));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
