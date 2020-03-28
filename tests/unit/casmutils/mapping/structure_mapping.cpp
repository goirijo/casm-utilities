// These are classes that structure_tools depends on
#include "../../../autotools.hh"
#include <gtest/gtest.h>

// This file tests the functions in:
#include <casmutils/mapping/structure_mapping.hpp>

namespace cu=casmutils;

class StructureMapTest : public testing::Test
{
protected:
    using Structure = casmutils::xtal::Structure;
    void SetUp() override
    {
        // Paths to testing poscars
        cu::fs::path primitive_fcc_path(cu::autotools::input_filesdir / "primitive_fcc_Ni.vasp");
        cu::fs::path primitive_bcc_path(cu::autotools::input_filesdir / "primitive_bcc_Ni.vasp");
        cu::fs::path partial_bain_path(cu::autotools::input_filesdir / "partial_bain_Ni.vasp");
        cu::fs::path perfect_bain_path(cu::autotools::input_filesdir / "perfectly_bained_Ni.vasp");
        cu::fs::path displaced_path(cu::autotools::input_filesdir / "displaced_fcc_Ni.vasp");
        // load the poscars
        primitive_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(primitive_fcc_path));
        primitive_bcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(primitive_bcc_path));
        partial_bain_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(partial_bain_path));
        perfect_bain_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(perfect_bain_path));
        displaced_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(displaced_path));
    }

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<cu::xtal::Structure> primitive_fcc_Ni_ptr;
    std::unique_ptr<cu::xtal::Structure> primitive_bcc_Ni_ptr;
    std::unique_ptr<cu::xtal::Structure> partial_bain_Ni_ptr;
    std::unique_ptr<cu::xtal::Structure> perfect_bain_Ni_ptr;
    std::unique_ptr<cu::xtal::Structure> displaced_fcc_Ni_ptr;
};

TEST_F(StructureMapTest, StructureMap)
{
    // map with fcc as reference and bcc,fully bained fcc and partially bained fcc as test structures
    cu::mapping::MappingReport full_bain =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *primitive_bcc_Ni_ptr)[0];
    cu::mapping::MappingReport perfect_bain =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *perfect_bain_Ni_ptr)[0];
    cu::mapping::MappingReport partial_bain =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *partial_bain_Ni_ptr)[0];

    auto [full_bain_lattice_score, full_bain_basis_score] = cu::mapping::structure_score(full_bain);
    auto [partial_bain_lattice_score, partial_bain_basis_score] = cu::mapping::structure_score(partial_bain);
    auto [perfect_bain_lattice_score, perfect_bain_basis_score] = cu::mapping::structure_score(perfect_bain);
    // lattice score should be finite and identical for bcc and fully bained no matter the volume
    EXPECT_TRUE(std::abs(full_bain_lattice_score - perfect_bain_lattice_score) < 1e-10);
    // partially bained fcc should have a lower score than bcc
    EXPECT_TRUE(full_bain_lattice_score > partial_bain_lattice_score);
    // basis scores should be 0
    EXPECT_TRUE(std::abs(full_bain_basis_score - partial_bain_basis_score) < 1e-10);
    EXPECT_TRUE(std::abs(perfect_bain_basis_score - partial_bain_basis_score) < 1e-10);
    EXPECT_TRUE(std::abs(perfect_bain_basis_score) < 1e-10);

    // map with fcc as reference and displaced fcc as test structure
    cu::mapping::MappingReport displaced =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *displaced_fcc_Ni_ptr)[0];
    auto [lattice_score, basis_score] = cu::mapping::structure_score(displaced);
    EXPECT_TRUE(std::abs(lattice_score) < 1e-10);
    EXPECT_TRUE(std::abs(basis_score - 0.08) < 1e-10);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
