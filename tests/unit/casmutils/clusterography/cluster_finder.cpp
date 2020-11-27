// These are classes that structure depends on
#include "../../../autotools.hh"
#include <casmutils/clusterography/cluster_finder.hpp>
#include <casmutils/xtal/structure.hpp>
#include <gtest/gtest.h>

class StructureTest : public testing::Test
{
protected:
    void SetUp() override {}

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<casmutils::xtal::Structure> cubic_Ni_struc_ptr;
};

TEST_F(StructureTest, ConstructandLatAccess)
{
    //  checks to see if Structure can be constructed
    // and has the lattice that was given associated with it.
    casmutils::fs::path pos_path(casmutils::autotools::input_filesdir / "conventional_fcc_Ni.vasp");
    casmutils::xtal::Structure pos = casmutils::xtal::Structure::from_poscar(pos_path);

    std::vector<double> max_length{3, 3, 8, 3, 3, 8};
    make_periodic_orbits(max_length, pos);

    EXPECT_TRUE(true);
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
