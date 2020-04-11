// These are classes that structure_tools depends on
#include "../../../autotools.hh"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

// This file tests the functions in:
#include <casmutils/mapping/structure_mapping.hpp>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cu = casmutils;

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

TEST_F(StructureMapTest, BainMappingScore)
{
    // map with fcc as reference and bcc,fully bained fcc and partially bained fcc as test structures
    cu::mapping::MappingReport full_bain_report = cu::mapping::map_structure(
        *primitive_fcc_Ni_ptr, *primitive_bcc_Ni_ptr)[0];
    cu::mapping::MappingReport perfect_bain_report = cu::mapping::map_structure(
        *primitive_fcc_Ni_ptr, *perfect_bain_Ni_ptr)[0];
    cu::mapping::MappingReport partial_bain_report = cu::mapping::map_structure(
        *primitive_fcc_Ni_ptr, *partial_bain_Ni_ptr)[0];

    auto [full_bain_lattice_score, full_bain_basis_score] =
        cu::mapping::structure_score(full_bain_report);
    auto [partial_bain_lattice_score, partial_bain_basis_score] =
        cu::mapping::structure_score(partial_bain_report);
    auto [perfect_bain_lattice_score, perfect_bain_basis_score] =
        cu::mapping::structure_score(perfect_bain_report);

    // lattice score should be finite and identical for bcc and fully bained no matter the volume
    EXPECT_TRUE(std::abs(full_bain_lattice_score - perfect_bain_lattice_score) < 1e-10);
    // partially bained fcc should have a lower score than bcc
    EXPECT_TRUE(full_bain_lattice_score > partial_bain_lattice_score);
    // basis scores should be 0
    EXPECT_TRUE(std::abs(full_bain_basis_score - partial_bain_basis_score) < 1e-10);
    EXPECT_TRUE(std::abs(perfect_bain_basis_score - partial_bain_basis_score) < 1e-10);
    EXPECT_TRUE(std::abs(perfect_bain_basis_score) < 1e-10);

}

TEST_F(StructureMapTest, DisplacementMappingScore)
{
    // map with fcc as reference and displaced fcc as test structure
    cu::mapping::MappingReport displacement_report = cu::mapping::map_structure(
        *primitive_fcc_Ni_ptr, *displaced_fcc_Ni_ptr)[0];
    auto [lattice_score, basis_score] =
        cu::mapping::structure_score(displacement_report);
    EXPECT_TRUE(std::abs(lattice_score) < 1e-10);
    EXPECT_TRUE(std::abs(basis_score - 0.08) < 1e-10);
}

//**********************************************************************************************

// TODO: Move this somewhere common because it's useful
std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

class GammaSurfaceMapTest : public testing::Test
{
protected:
    using Structure = casmutils::xtal::Structure;
    void SetUp() override
    {
        root = cu::autotools::input_filesdir / "Mg-mush/basal";
        /* root = cu::autotools::input_filesdir / "Li-mush/111"; */
        auto zero_cleave_dirs = this->directories_zero_cleave();

        auto cmp_by_coord = [](const cu::fs::path& lhs, const cu::fs::path& rhs) {
            return coordinate_from_directory(lhs) < coordinate_from_directory(rhs);
        };
        std::sort(zero_cleave_dirs.begin(), zero_cleave_dirs.end(), cmp_by_coord);

        for (const auto& cleave_dir : zero_cleave_dirs)
        {
            shifted_coordinates.emplace_back(coordinate_from_directory(cleave_dir));
            shifted_structures.emplace_back(cu::xtal::Structure::from_poscar(cleave_dir / "POSCAR"));
        }

        return;
    }

    std::vector<cu::xtal::Structure> shifted_structures;
    std::vector<std::pair<int, int>> shifted_coordinates;

    cu::fs::path root;

private:
    std::vector<cu::fs::path> directories_zero_cleave() const
    {
        std::vector<cu::fs::path> no_cleave_dirs;
        for (const auto& d : cu::fs::recursive_directory_iterator(this->root))
        {
            if (!d.is_directory())
            {
                continue;
            }

            cu::fs::path dir(d);
            if (dir.filename() != "cleave_0.000000")
            {
                continue;
            }

            no_cleave_dirs.emplace_back(std::move(dir));
        }

        return no_cleave_dirs;
    }

    static std::pair<int, int> coordinate_from_directory(const cu::fs::path cleave_dir)
    {
        int x, y;

        std::string shift_dirname = cleave_dir.parent_path().filename();
        std::vector<std::string> split_by_period = split(shift_dirname, '.');
        y = std::stoi(split_by_period.back());

        std::vector<std::string> split_by_underscore = split(split_by_period[0], '_');
        x = std::stoi(split_by_underscore.back());

        return std::make_pair(x, y);
    }
};

TEST_F(GammaSurfaceMapTest, MapAllShifts) {
    //Gamma surface of basal plane has mirror symmetry, so expect
    //half of the structures to be equivalent to the other half
    cu::mapping::MappingInput map_strategy;
    /* map_strategy.use_crystal_symmetry=true; */
    map_strategy.k_best_maps=0;
    map_strategy.min_cost=1e-5;

    cu::xtal::Structure reference=shifted_structures[0];

    for(int i=0; i<shifted_structures.size(); ++i)
    {
        std::cout<<i<<" "<<shifted_coordinates[i].first<<" "<<shifted_coordinates[i].second<<" ";



        cu::mapping::StructureMapper_f map_to_ith_struc(shifted_structures[i],map_strategy);

        for(int j=0; j<shifted_structures.size(); ++j)
        {
            auto reports=map_to_ith_struc(shifted_structures[j]);
            if(reports.size())
            {
                std::cout<<j<<" ";
                /* std::cout<<"\n"<<reports[0].cost<<" "<<reports[0].lattice_cost<<" "<<reports[0].cost<<"\n"; */
            }
            /* auto [l,b]=cu::mapping::structure_score(reports[0]); */
            /* std::cout<<l<<","<<b<<std::endl; */
        }
        std::cout<<std::endl;
    }

    //For every map, make sure that all structures in the map agree that they're the same
    //
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
