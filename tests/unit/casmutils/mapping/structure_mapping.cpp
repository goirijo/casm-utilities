// These are classes that structure_tools depends on
#include "../../../autotools.hh"
#include <algorithm>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

// This file tests the functions in:
#include <casmutils/mapping/structure_mapping.hpp>
#include <limits>
#include <math.h>
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
    cu::mapping::MappingReport full_bain_report =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *primitive_bcc_Ni_ptr)[0];
    cu::mapping::MappingReport perfect_bain_report =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *perfect_bain_Ni_ptr)[0];
    cu::mapping::MappingReport partial_bain_report =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *partial_bain_Ni_ptr)[0];

    auto [full_bain_lattice_score, full_bain_basis_score] = cu::mapping::structure_score(full_bain_report);
    auto [partial_bain_lattice_score, partial_bain_basis_score] = cu::mapping::structure_score(partial_bain_report);
    auto [perfect_bain_lattice_score, perfect_bain_basis_score] = cu::mapping::structure_score(perfect_bain_report);

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
    cu::mapping::MappingReport displacement_report =
        cu::mapping::map_structure(*primitive_fcc_Ni_ptr, *displaced_fcc_Ni_ptr)[0];
    auto [lattice_score, basis_score] = cu::mapping::structure_score(displacement_report);
    EXPECT_TRUE(std::abs(lattice_score) < 1e-10);
    EXPECT_TRUE(std::abs(basis_score - 0.08) < 1e-10);
}

class SymmetryPreservingMappingTest : public testing::Test
{
protected:
    using Structure = casmutils::xtal::Structure;
    void SetUp() override
    {

        cu::fs::path tall_hcp_path(cu::autotools::input_filesdir / "tall_hcp.vasp");
        cu::fs::path squished_hcp_path(cu::autotools::input_filesdir / "squished_hcp.vasp");

        cu::fs::path conventional_fcc_path(cu::autotools::input_filesdir / "conventional_fcc_Ni.vasp");
        cu::fs::path partial_bain_path(cu::autotools::input_filesdir / "partial_bain_Ni.vasp");

        tall_hcp_ptr = std::make_unique<Structure>(Structure::from_poscar(tall_hcp_path));
        squished_hcp_ptr = std::make_unique<Structure>(Structure::from_poscar(squished_hcp_path));

        conventional_fcc_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(conventional_fcc_path));
        partial_bain_Ni_ptr = std::make_unique<Structure>(Structure::from_poscar(partial_bain_path));
    }

    std::unique_ptr<Structure> tall_hcp_ptr;
    std::unique_ptr<Structure> squished_hcp_ptr;
    std::unique_ptr<Structure> conventional_fcc_Ni_ptr;
    std::unique_ptr<Structure> partial_bain_Ni_ptr;
    double tol = 1e-5;
};

TEST_F(SymmetryPreservingMappingTest, PreservingTest)
{
    cu::mapping::MappingReport full_report = cu::mapping::map_structure(*tall_hcp_ptr, *squished_hcp_ptr)[0];
    auto hcp_group = cu::xtal::make_factor_group(*tall_hcp_ptr, tol);
    // construct the corresponding permutation representation
    cu::sym::PermRep no_swap = {0, 1};
    cu::sym::PermRep swap = {1, 0};
    std::vector<cu::sym::PermRep> perm_group;
    for (const auto op : hcp_group)
    {
        bool op_swaps = false;
        if (op_swaps)
        {
            perm_group.push_back(swap);
        }
        else
        {
            perm_group.push_back(no_swap);
        }
    }
    cu::mapping::MappingReport adjusted_report =
        cu::mapping::symmetry_preserving_mapping_report(full_report, hcp_group, perm_group);
    // because only difference is c/a ratio the mapping report should be entirely symmetry preserving
    EXPECT_TRUE(cu::is_equal(full_report.displacement, adjusted_report.displacement, 1e-5));
    EXPECT_TRUE(cu::is_equal(full_report.stretch, adjusted_report.stretch, 1e-5));
}

TEST_F(SymmetryPreservingMappingTest, NotPreservingTest)
{
    cu::mapping::MappingReport full_report =
        cu::mapping::map_structure(*conventional_fcc_Ni_ptr, *partial_bain_Ni_ptr)[0];
    auto fcc_group = cu::xtal::make_factor_group(*conventional_fcc_Ni_ptr, tol);
    std::cout << "DEBUGGING: full_report.stretch" << full_report.stretch << std::endl;
    // construct the corresponding permutation representation
    cu::sym::PermRep no_swap = {0, 1, 2, 3};
    cu::sym::PermRep swap = {3, 2, 1, 0};
    std::vector<cu::sym::PermRep> perm_group;
    for (const auto op : fcc_group)
    {
        bool op_swaps = false;
        if (op_swaps)
        {
            perm_group.push_back(swap);
        }
        else
        {
            perm_group.push_back(no_swap);
        }
    }
    cu::mapping::MappingReport adjusted_report =
        cu::mapping::symmetry_preserving_mapping_report(full_report, fcc_group, perm_group);
    double volume_factor = ((conventional_fcc_Ni_ptr->lattice().volume()) / partial_bain_Ni_ptr->lattice().volume());
    double scale_factor = (volume_factor - 1.0) / 3 + 1;
    // because only difference is partial bain ratio the mapping report should be entirely symmetry breaking
    EXPECT_TRUE(cu::is_equal(full_report.displacement, adjusted_report.displacement, 1e-5));
    Eigen::Matrix3d ident = Eigen::Matrix3d::Identity() * scale_factor;
    EXPECT_TRUE(cu::is_equal(ident, adjusted_report.stretch, 1e-5));
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

class MgGammaSurfaceMapTest : public testing::Test
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

        a_shifts = b_shifts = 0;
        for (const auto& cleave_dir : zero_cleave_dirs)
        {
            shifted_coordinates.emplace_back(coordinate_from_directory(cleave_dir));
            shifted_structures.emplace_back(cu::xtal::Structure::from_poscar(cleave_dir / "POSCAR"));

            a_shifts = std::max(a_shifts, shifted_coordinates.back().first + 1);
            b_shifts = std::max(b_shifts, shifted_coordinates.back().second + 1);
        }
        assert(a_shifts == b_shifts);

        this->fill_equivalence_grid();
        return;
    }

    std::vector<cu::xtal::Structure> shifted_structures;
    std::vector<std::pair<int, int>> shifted_coordinates;

    int a_shifts, b_shifts;

    cu::fs::path root;

    // Store the indexes of the equivalent structures on the shift grid
    std::vector<std::vector<std::vector<int>>> equivalence_grid;

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

    void fill_equivalence_grid()
    {
        cu::mapping::MappingInput map_strategy;
        map_strategy.use_crystal_symmetry = true;
        map_strategy.k_best_maps = 0;
        map_strategy.min_cost = 1e-10;

        equivalence_grid =
            std::vector<std::vector<std::vector<int>>>(a_shifts, std::vector<std::vector<int>>(b_shifts));

        for (int i = 0; i < shifted_structures.size(); ++i)
        {
            cu::mapping::StructureMapper_f map_to_ith_struc(shifted_structures[i], map_strategy);

            auto [a, b] = shifted_coordinates[i];
            for (int j = 0; j < shifted_structures.size(); ++j)
            {
                auto reports = map_to_ith_struc(shifted_structures[j]);
                if (reports.size())
                {
                    equivalence_grid[a][b].push_back(j);
                }
            }
        }
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

TEST_F(MgGammaSurfaceMapTest, MirrorMappingSymmetry)
{
    // Gamma surface of basal plane has mirror symmetry, so expect
    // half of the structures to be equivalent to the other half

    // Expect at least mirror symmetry, since for HCP, shifting along a is the same as
    // shifting along b, provided the increments are the same, AND the vectors aren't obtuse
    for (int i = 0; i < equivalence_grid.size(); ++i)
    {
        for (int j = 0; j < equivalence_grid[i].size(); ++j)
        {
            EXPECT_EQ(equivalence_grid[i][j], equivalence_grid[j][i]);
        }
    }
}

TEST_F(MgGammaSurfaceMapTest, FactorGroupFactorMapping)
{
    // The number of equivalent structures should be a factor of the factor group size (24)
    for (int i = 0; i < equivalence_grid.size(); ++i)
    {
        for (int j = 0; j < equivalence_grid[i].size(); ++j)
        {
            int maps = equivalence_grid[i][j].size();
            EXPECT_TRUE(maps == 1 || maps == 2 || maps == 3 || maps == 4 || maps == 6);
        }
    }
}

TEST_F(MgGammaSurfaceMapTest, AtLeastSomeSymmetry)
{
    // Make sure it's not just all structures mapping only onto themselves
    bool mapped_6 = false;
    for (int i = 0; i < equivalence_grid.size(); ++i)
    {
        for (int j = 0; j < equivalence_grid[i].size(); ++j)
        {
            if (equivalence_grid[i][j].size() == 6)
            {
                mapped_6 = true;
            }
        }
    }

    EXPECT_TRUE(mapped_6);
}

TEST_F(MgGammaSurfaceMapTest, MappingResultsSize)
{
    cu::xtal::Structure hcp_triple = shifted_structures[0];
    auto factor_group = cu::xtal::make_factor_group(hcp_triple, 1e-5);
    // If you don't use symmetry in the mapper, expect to get as many
    // mappings as there are factor group operations
    cu::mapping::MappingInput map_strategy;
    map_strategy.k_best_maps = 0;
    map_strategy.min_cost = 1e-10;

    cu::mapping::StructureMapper_f blind_mapper(hcp_triple, map_strategy);
    EXPECT_EQ(blind_mapper(hcp_triple).size(), factor_group.size());

    cu::mapping::StructureMapper_f sym_aware_mapper(hcp_triple, map_strategy, factor_group);
    EXPECT_EQ(sym_aware_mapper(hcp_triple).size(), 3);
}

TEST_F(MgGammaSurfaceMapTest, SelfMapResultsSize)
{
    const auto& hcp_super = this->shifted_structures[0];

    cu::mapping::MappingInput map_strategy;
    map_strategy.use_crystal_symmetry = true;
    map_strategy.k_best_maps = 0;
    map_strategy.min_cost = 1e-10;

    cu::mapping::StructureMapper_f map_to_hcp_super_with_sym(hcp_super, map_strategy);
    // When using crystal symmetry, the mapper should only find as many mappings as primitives
    // fit in the structure
    EXPECT_EQ(3, map_to_hcp_super_with_sym(hcp_super).size());

    map_strategy.use_crystal_symmetry = false;
    cu::mapping::StructureMapper_f map_to_hcp_super_without_sym(hcp_super, map_strategy);
    // When not using crystal symmetry, the mapper should find as many mappings as
    // there are factor group operations (3*24, i.e. number of primitives * hcp factor group size)
    EXPECT_EQ(72, map_to_hcp_super_without_sym(hcp_super).size());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
