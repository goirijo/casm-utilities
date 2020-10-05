#include "../../../autotools.hh"
#include "casmutils/mush/slab.hpp"
#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>

#include "../../../autotools.hh"
#include <casmutils/definitions.hpp>
#include <casmutils/mush/twist.hpp>
#include <casmutils/xtal/frankenstein.hpp>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace cu = casmutils;
using Lattice = cu::xtal::Lattice;

class TwistTest : public testing::Test
{
protected:
    void SetUp()
    {
        Eigen::Matrix3d col_lat_mat;
        col_lat_mat << 3, 0, 0, 0, 4, 0, 0, 0, 8;
        Lattice ortho_lat(col_lat_mat.col(0), col_lat_mat.col(1), col_lat_mat.col(2));

        for (int i : {-1, 0, 2})
        {
            for (int j : {-2, 1, 3})
            {
                for (int k : {1, -4, 2})
                {
                    Eigen::Vector3i millers(i, j, k);
                    sliced_lattices.emplace_back(cu::xtal::slice_along_plane(ortho_lat, millers));
                }
            }
        }
        return;
    }

    std::vector<Lattice> sliced_lattices;
    double tol = 1e-10;

    // TODO: Move this somewhere useful
    static double degrees_between_vectors(const Eigen::Vector3d& a1, const Eigen::Vector3d& a2)
    {
        return 2 * std::atan2((a2.norm() * a1 - a1.norm() * a2).norm(), (a2.norm() * a1 + a1.norm() * a2).norm()) *
               180 / M_PI;
    }

    static double
    signed_degrees_between_vectors(const Eigen::Vector3d& a1, const Eigen::Vector3d& a2, const Eigen::Vector3d& normal)
    {
        int sign = 1;
        if (normal.dot(a1.cross(a2)) < 0)
        {
            sign = -1;
        }
        return sign * degrees_between_vectors(a1, a2);
    }
};

TEST_F(TwistTest, TwistMatrixIsUnitary)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double degrees : {-5, 0, 30})
        {

            Eigen::Matrix3d R = cu::mush::make_twist_rotation_matrix(start_lat, degrees);
            EXPECT_TRUE(cu::almost_equal(R.determinant(), 1.0, tol));

            // Unitary if R.T*R is identity
            Eigen::Matrix3d RtransposeR = R.transpose() * R;
            Eigen::Matrix3d I(Eigen::Matrix3d::Identity());

            EXPECT_TRUE(almost_equal(RtransposeR, I, tol));
        }
    }
}

TEST_F(TwistTest, TwistedLatticeAngles)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        Eigen::Vector3d normal = start_lat.a().cross(start_lat.b());
        for (double angle : {-33.0, -15.0, -1.0, -0.1, 0.0, 0.00000001, 0.2, 30.0, 90.0, 180.0})
        {
            Lattice rotated_lat = cu::mush::make_twisted_lattice(start_lat, angle);

            EXPECT_TRUE(cu::almost_equal(normal, rotated_lat.a().cross(rotated_lat.b()), tol));
            EXPECT_TRUE(cu::almost_equal(start_lat.c().dot(normal), rotated_lat.c().dot(normal), tol));

            double a_deg = signed_degrees_between_vectors(start_lat.a(), rotated_lat.a(), normal);
            double b_deg = signed_degrees_between_vectors(start_lat.b(), rotated_lat.b(), normal);

            EXPECT_TRUE(cu::almost_equal(a_deg, b_deg, tol)) << a_deg << " vs " << b_deg;
            EXPECT_TRUE(cu::almost_equal(a_deg, angle, tol)) << a_deg << " vs " << angle;
        }
    }
}

TEST_F(TwistTest, ObviousTwistMatrix)
{
    Eigen::Matrix3d ortho_mat;
    ortho_mat << 4, 0, 0, 0, 6, 0, 0, 0, 9;

    cu::xtal::Lattice ortho_lat(ortho_mat);

    Eigen::Matrix3d rotate_0 = cu::mush::make_twist_rotation_matrix(ortho_lat, 0);
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());
    EXPECT_TRUE(almost_equal(rotate_0, I, tol));

    Eigen::Matrix3d rotate_90 = cu::mush::make_twist_rotation_matrix(ortho_lat, 90);
    Eigen::Matrix3d exptected_rotate_90;
    exptected_rotate_90 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    EXPECT_TRUE(almost_equal(rotate_90, exptected_rotate_90, 1e-10));
}

TEST_F(TwistTest, LessObviousTwistMatrix)
{
    Eigen::Matrix3d ortho_mat;
    ortho_mat << 0, 0, 9, 4, 0, 0, 0, 6, 0;

    cu::xtal::Lattice ortho_lat(ortho_mat);

    Eigen::Matrix3d rotate_0 = cu::mush::make_twist_rotation_matrix(ortho_lat, 0);
    Eigen::Matrix3d I(Eigen::Matrix3d::Identity());
    EXPECT_TRUE(almost_equal(rotate_0, I, 1e-10));

    Eigen::Matrix3d rotate_90 = cu::mush::make_twist_rotation_matrix(ortho_lat, 90);
    Eigen::Matrix3d exptected_rotate_90;
    exptected_rotate_90 << 1, 0, 0, 0, 0, -1, 0, 1, 0;

    EXPECT_TRUE(almost_equal(rotate_90, exptected_rotate_90, tol));
}

TEST_F(TwistTest, AlginLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        cu::xtal::Lattice aligned_lat = cu::mush::make_aligned_lattice(start_lat);
        EXPECT_TRUE(cu::almost_equal(aligned_lat.a()(1), 0.0, tol));
        EXPECT_TRUE(cu::almost_equal(aligned_lat.a()(2), 0.0, tol));
        EXPECT_TRUE(cu::almost_equal(aligned_lat.b()(2), 0.0, tol));

        EXPECT_TRUE(cu::almost_equal(start_lat.a().norm(), aligned_lat.a().norm(), tol));
        EXPECT_TRUE(cu::almost_equal(start_lat.b().norm(), aligned_lat.b().norm(), tol));
        EXPECT_TRUE(cu::almost_equal(start_lat.c().norm(), aligned_lat.c().norm(), tol));
    }
}

TEST_F(TwistTest, ApproximantMoireLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE = cu::mush::MoireApproximant::LATTICE;
            cu::mush::MoireLattice moire(start_lat, angle);
            cu::mush::MoireApproximant approx_moire(
                moire.aligned_moire_lattice, moire.aligned_lattice, moire.rotated_lattice);
            const auto& moire_lat = approx_moire.approximate_moire_lattice;
            const auto& aligned_lat = approx_moire.approximate_lattices.at(LATTICE::ALIGNED);
            const auto& rot_lat = approx_moire.approximate_lattices.at(LATTICE::ROTATED);

            Eigen::Matrix3d aligned_to_moire_transform =
                aligned_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff =
                aligned_to_moire_transform - aligned_to_moire_transform.unaryExpr([](double x) {
                    return std::round(x);
                });

            Eigen::Matrix3d rot_to_moire_transform =
                rot_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff =
                rot_to_moire_transform - rot_to_moire_transform.unaryExpr([](double x) {
                    return std::round(x);
                });

            EXPECT_TRUE(almost_zero(aligned_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
            EXPECT_TRUE(almost_zero(rot_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
        }
    }
}

TEST_F(TwistTest, MoireApproximantStrainMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE = cu::mush::MoireApproximant::LATTICE;
            cu::mush::MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            cu::mush::MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            for (LATTICE lat : {LATTICE::ALIGNED, LATTICE::ROTATED})
            {
                const Eigen::Matrix3d& L = moire.real(lat).column_vector_matrix();
                const Eigen::Matrix3d& L_prime = approx_moire.approximate_lattices.at(lat).column_vector_matrix();

                Eigen::Matrix3d F = L_prime * L.inverse();
                EXPECT_TRUE(almost_equal(F, approx_moire.approximation_deformations.at(lat)));
            }
        }
    }
}

TEST_F(TwistTest, MoireApproximantVectorMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE = cu::mush::MoireApproximant::LATTICE;
            cu::mush::MoireLattice moire(start_lat, angle);
            const Lattice& moire_lat = moire.aligned_moire_lattice;
            const Lattice& aligned_lat = moire.aligned_lattice;
            const Lattice& rotated_lat = moire.rotated_lattice;

            cu::mush::MoireApproximant approx_moire(moire_lat, aligned_lat, rotated_lat);

            Lattice aligned_super = cu::xtal::make_superlattice(
                approx_moire.approximate_lattices.at(LATTICE::ALIGNED),
                approx_moire.approximate_moire_integer_transformations.at(LATTICE::ALIGNED).cast<int>());

            Lattice rot_super = cu::xtal::make_superlattice(
                approx_moire.approximate_lattices.at(LATTICE::ROTATED),
                approx_moire.approximate_moire_integer_transformations.at(LATTICE::ROTATED).cast<int>());

            EXPECT_TRUE(almost_equal(aligned_super.a(), rot_super.a()));
            EXPECT_TRUE(almost_equal(aligned_super.b(), rot_super.b()));
        }
    }
}

TEST_F(TwistTest, MoireApproximatorVectorMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            cu::mush::MoireApproximator approx_moire(start_lat, angle);
            using ZONE = cu::mush::MoireApproximator::ZONE;
            using LATTICE = cu::mush::MoireApproximator::LATTICE;

            const auto approx_moire_AA = approx_moire.best_smallest(ZONE::ALIGNED, LATTICE::ALIGNED);
            Lattice aligned_super = cu::xtal::make_superlattice(
                approx_moire_AA.approximate_tiling_unit, approx_moire_AA.tiling_unit_supercell_matrix.cast<int>());

            const auto approx_moire_AR = approx_moire.best_smallest(ZONE::ALIGNED, LATTICE::ROTATED);
            Lattice rot_super = cu::xtal::make_superlattice(approx_moire_AR.approximate_tiling_unit,
                                                            approx_moire_AR.tiling_unit_supercell_matrix.cast<int>());

            EXPECT_TRUE(almost_equal(aligned_super.a(), rot_super.a()));
            EXPECT_TRUE(almost_equal(aligned_super.b(), rot_super.b()));
        }
    }
}

TEST_F(TwistTest, MoireApproximantAlignedRotatedMatch)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            // TODO: Make sure that the aligned and rotated lattices
            //(both perfect and approximated) are exclusively related
            // by a rotation
        }
    }
}

TEST_F(TwistTest, ApproximantPrismaticMoireDeformation)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            using LATTICE = cu::mush::MoireLattice::LATTICE;
            cu::mush::MoireLattice moire(start_lat, angle);
            cu::mush::MoireApproximant approx_moire(
                moire.aligned_moire_lattice, moire.aligned_lattice, moire.rotated_lattice);
            const auto& moire_lat = approx_moire.approximate_moire_lattice;
            const auto& aligned_lat = approx_moire.approximate_lattices.at(LATTICE::ALIGNED);
            const auto& rot_lat = approx_moire.approximate_lattices.at(LATTICE::ROTATED);

            Eigen::Matrix3d aligned_to_moire_transform =
                aligned_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d aligned_to_moire_transform_diff =
                aligned_to_moire_transform - aligned_to_moire_transform.unaryExpr([](double x) {
                    return std::round(x);
                });

            Eigen::Matrix3d rot_to_moire_transform =
                rot_lat.column_vector_matrix().inverse() * moire_lat.column_vector_matrix();
            Eigen::Matrix3d rot_to_moire_transform_diff =
                rot_to_moire_transform - rot_to_moire_transform.unaryExpr([](double x) {
                    return std::round(x);
                });

            EXPECT_TRUE(almost_zero(aligned_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
            EXPECT_TRUE(almost_zero(rot_to_moire_transform_diff.block<2, 2>(0, 0), 1e-8));
        }
    }
}

TEST_F(TwistTest, PrismaticLattice)
{
    for (const Lattice& start_lat : sliced_lattices)
    {
        for (double angle : {-3.00, -2.0, -1.0, -0.50, 0.50, 2.00, 5.00, 22.0})
        {
            Lattice prismatic_lat = cu::mush::make_prismatic_lattice(start_lat);
            EXPECT_EQ(start_lat.a(), prismatic_lat.a());
            EXPECT_EQ(start_lat.b(), prismatic_lat.b());

            EXPECT_TRUE(cu::almost_equal(prismatic_lat.a().dot(prismatic_lat.c()), 0.0));
            EXPECT_TRUE(cu::almost_equal(prismatic_lat.b().dot(prismatic_lat.c()), 0.0));
        }
    }
}

TEST(HexagonalTwistTest, BrillouinOverlapMoiresRelatedByIntegerTransform)
{
    int passed = 0;
    int missed = 0;
    Eigen::Matrix3d row_lat_mat;
    row_lat_mat << 2.4684159756, 0.0000000000, 0.0000000000, -1.2342079878, 2.1377109420, 0.0000000000, 0.0000000000,
        0.0000000000, 9.9990577698;

    Lattice triangular(row_lat_mat.row(0), row_lat_mat.row(1), row_lat_mat.row(2));

    for (double degrees = 0.5; degrees < 361; degrees += 1.0)
    {
        bool full_brillouin_overlap = true;

        cu::mush::MoireLattice moire(triangular, degrees);
        for (auto lat : {cu::mush::MoireLattice::LATTICE::ALIGNED, cu::mush::MoireLattice::LATTICE::ROTATED})
        {
            for (int i = 0; i < 2; ++i)
            {
                if (!moire.is_within_brillouin_zone_overlap[lat][i])
                {
                    full_brillouin_overlap = false;
                }
            }
        }

        EXPECT_TRUE(full_brillouin_overlap);

        // TODO: Expose transformation matrix stuff (Superlattice?) in cu
        const auto& M = moire.aligned_moire_lattice;
        const auto& Mt = moire.rotated_moire_lattice;
        Eigen::Matrix3d m_to_m = M.column_vector_matrix().inverse() * Mt.column_vector_matrix();

        Eigen::Matrix3i im_to_m = cu::iround(m_to_m);
        Eigen::Matrix3d error = m_to_m - im_to_m.cast<double>();

        EXPECT_TRUE(almost_equal(error, Eigen::Matrix3d::Zero()));

        ++passed;
    }
}

class GrapheneTwistTest : public testing::Test
{
public:
    using LAT = cu::mush::MoireApproximator::LATTICE;
    using ZONE = cu::mush::MoireApproximator::ZONE;
    using Structure = cu::xtal::Structure;

protected:
    void SetUp()
    {
        graphene_ptr.reset(new Structure(Structure::from_poscar(cu::autotools::input_filesdir / "graphene.vasp")));

        magic_angles.push_back(2.449977276616);
        magic_angles.push_back(2.645908381192);
        magic_angles.push_back(2.875894633632);
        magic_angles.push_back(3.149657426389);
        magic_angles.push_back(3.481006089466);
        magic_angles.push_back(3.890238169007);
        magic_angles.push_back(4.408455007944);
        magic_angles.push_back(5.085847808123);
        magic_angles.push_back(5.509040771928);
        magic_angles.push_back(6.008983197766);
        magic_angles.push_back(6.608610360312);
        magic_angles.push_back(7.340993016630);
        magic_angles.push_back(7.926469934899);
        magic_angles.push_back(8.255620609025);
        magic_angles.push_back(8.613238191002);
        magic_angles.push_back(9.430007907896);
        magic_angles.push_back(10.41743820571);
        magic_angles.push_back(10.99273308930);
        magic_angles.push_back(11.63505128888);
        magic_angles.push_back(11.98510029238);
        magic_angles.push_back(13.17355110725);
        magic_angles.push_back(14.30767637471);
        magic_angles.push_back(14.62222138876);
        magic_angles.push_back(15.17817893794);
        magic_angles.push_back(15.65414359251);
        magic_angles.push_back(16.42642140347);
        magic_angles.push_back(17.27824434988);
        magic_angles.push_back(17.89655112925);
        magic_angles.push_back(18.73399783335);
        magic_angles.push_back(19.27480690356);
        magic_angles.push_back(19.65285963166);
        magic_angles.push_back(21.78678929826);
        magic_angles.push_back(24.01660222864);
        magic_angles.push_back(24.43269767945);
        magic_angles.push_back(25.03965959447);
        magic_angles.push_back(26.00782388564);
        magic_angles.push_back(26.74565161623);
        magic_angles.push_back(27.79577249602);
        magic_angles.push_back(28.78320279384);
        magic_angles.push_back(29.40931139719);
        magic_angles.push_back(30.15827583531);
        magic_angles.push_back(30.59068860280);
        magic_angles.push_back(32.20422750397);
        magic_angles.push_back(33.99217611435);
        magic_angles.push_back(34.53894442320);
        magic_angles.push_back(35.56730232054);
        magic_angles.push_back(36.51693800821);
        magic_angles.push_back(38.21321070173);
        magic_angles.push_back(39.68334042751);
        magic_angles.push_back(40.34714036833);
        magic_angles.push_back(40.96932360328);
        magic_angles.push_back(42.10344887074);
        magic_angles.push_back(43.57357859652);
        magic_angles.push_back(44.82182106205);
        magic_angles.push_back(45.89464594547);
        magic_angles.push_back(46.82644889274);
        magic_angles.push_back(48.36494871111);
        magic_angles.push_back(49.58256179428);
        magic_angles.push_back(50.56999209210);
        magic_angles.push_back(51.38676180899);
        magic_angles.push_back(52.07353006510);
        magic_angles.push_back(52.65900698336);
        magic_angles.push_back(53.16403781617);
    }

    std::unique_ptr<cu::xtal::Structure> graphene_ptr;
    // Angles that give coindicent moire superlattices, but may require supercells.
    std::vector<double> magic_angles;
};

TEST_F(GrapheneTwistTest, MoireScel15DegreeTwist)
{
    /* auto equivalent=[](const cu::xtal::Lattice& L1, const cu::xtal::Lattice& L2) */
    /* { */
    /*     auto [T,E]=cu::mush::approximate_integer_transformation(L1,L2); */
    /*     return almost_zero(E) && CASM::is_unimodular(T,1e-8); */
    /* }; */

    cu::xtal::LatticeIsEquivalent_f equivalent(1e-8);

    double angle = 15.178178937949879; // This angle gives a coincident sqrt(3) sqrt(3) moire superlattice
    Eigen::Matrix3i sqrt3_transfmat;
    sqrt3_transfmat << 2, 1, 0, 1, 2, 0, 0, 0, 1;

    cu::mush::MoireApproximator mini_graph_moire(graphene_ptr->lattice(), angle);
    const auto& mini_moire_unit = mini_graph_moire.true_moire(ZONE::ALIGNED);
    const auto mini_moire_best_smallest=mini_graph_moire.best_smallest(ZONE::ALIGNED,LAT::ALIGNED,1e-10);
    mini_graph_moire.expand(100);
    const auto better_moire_best_smallest=mini_graph_moire.best_smallest(ZONE::ALIGNED,LAT::ALIGNED,1e-10);

    cu::xtal::Lattice sqrt3_super_moire = cu::xtal::make_superlattice(mini_moire_unit, sqrt3_transfmat);

    EXPECT_TRUE(equivalent(mini_moire_best_smallest.true_moire,mini_moire_unit));
    EXPECT_TRUE(equivalent(better_moire_best_smallest.approximate_moire, sqrt3_super_moire));
}

TEST_F(GrapheneTwistTest, MoireScelMagicDegreeTwist)
{
    const auto I = Eigen::Matrix3d::Identity();
    for (double angle : magic_angles)
    {
        cu::mush::MoireApproximator graph_moire(graphene_ptr->lattice(), angle);
        graph_moire.expand(1000);
        for (auto Z : {ZONE::ALIGNED, ZONE::ROTATED})
        {
            for (auto L : {LAT::ALIGNED, LAT::ROTATED})
            {
                const auto report = graph_moire.best_smallest(Z, L, 1e-8);
                const auto& F = report.approximation_deformation;
                EXPECT_TRUE(almost_zero(F - I));
            }
        }
    }
}

void spit_rotation_error_data(const cu::xtal::Structure& tile, const std::vector<double>& magic_angles, int max_tiles)
{
    using ZONE = cu::mush::MoireApproximator::ZONE;
    using LAT = cu::mush::MoireApproximator::LATTICE;

    std::ofstream incdeg("size" + std::to_string(max_tiles) + ".txt");
    for (double degrees = -60.25; degrees < 61; degrees += 0.5)
    {
        cu::mush::MoireApproximator graph_moire(tile.lattice(), degrees, max_tiles);
        const auto aligned_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED);
        const auto rotated_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED);

        incdeg << degrees << "    ";
        for (auto report : {aligned_report, rotated_report})
        {
            cu::mush::DeformationReport dr(report.approximation_deformation);
            for (auto mat : {dr.rotation, dr.strain})
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        incdeg << mat(i, j) << "    ";
                    }
                }
            }
        }
        incdeg << "\n";
    }
    incdeg.close();

    std::ofstream magicdeg("magic" + std::to_string(max_tiles) + ".txt");
    auto double_magic = magic_angles;

    for (double d : magic_angles)
    {
        double_magic.push_back(-d);
    }

    std::ofstream magicdet("magic" + std::to_string(max_tiles) + "det.txt");
    for (double degrees : double_magic)
    {
        cu::mush::MoireApproximator graph_moire(tile.lattice(), degrees, max_tiles);
        const auto aligned_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED);
        const auto rotated_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED);

        magicdet << degrees << "    " << aligned_report.tiling_unit_supercell_matrix.determinant() << "    "
                 << rotated_report.tiling_unit_supercell_matrix.determinant() << "    "
                 << aligned_report.true_moire_supercell_matrix.determinant() << "    "
                 << rotated_report.true_moire_supercell_matrix.determinant() << "\n";

        magicdeg << degrees << "    ";
        for (auto report : {aligned_report, rotated_report})
        {
            cu::mush::DeformationReport dr(report.approximation_deformation);
            for (auto mat : {dr.rotation, dr.strain})
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        magicdeg << mat(i, j) << "    ";
                    }
                }
            }
        }
        magicdeg << "\n";
    }
    magicdeg.close();

    return;
}

void spit_determinants(const cu::xtal::Lattice& tile, const std::vector<double>& angles, int max_tiles)
{
    using ZONE = cu::mush::MoireApproximator::ZONE;
    using LAT = cu::mush::MoireApproximator::LATTICE;

    std::ofstream magicdet("magic" + std::to_string(max_tiles) + "det.txt");
    for (double degrees : angles)
    {
        cu::mush::MoireApproximator graph_moire(tile, degrees, max_tiles);
        const auto aligned_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED);
        const auto rotated_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED);

        magicdet << degrees << "    " << aligned_report.tiling_unit_supercell_matrix.determinant() << "    "
                 << rotated_report.tiling_unit_supercell_matrix.determinant() << "    "
                 << aligned_report.true_moire_supercell_matrix.determinant() << "    "
                 << rotated_report.true_moire_supercell_matrix.determinant() << "\n";

    }
    magicdet.close();
}

void spit_internal_error_data(const cu::xtal::Lattice& tile, const std::vector<double>& magic_angles, int max_tiles)
{
    using ZONE = cu::mush::MoireApproximator::ZONE;
    using LAT = cu::mush::MoireApproximator::LATTICE;

    std::ofstream incdeg("internal_size" + std::to_string(max_tiles) + ".txt");
    for (double degrees = -60.25; degrees < 61; degrees += 0.5)
    {
        cu::mush::MoireApproximator graph_moire(tile, degrees, max_tiles);
        const auto aligned_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED);
        const auto rotated_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED);

        incdeg << degrees << "    ";
        for (auto report : {aligned_report, rotated_report})
        {
            cu::mush::DeformationReport dr(report.approximation_deformation);
            incdeg << dr.rotation_angle << "    " << dr.dilation_strain << "    " << dr.deviatoric_strain << "    ";
        }
        incdeg << "\n";
    }
    incdeg.close();

    std::ofstream magicdeg("internal_magic" + std::to_string(max_tiles) + ".txt");
    auto double_magic = magic_angles;

    for (double d : magic_angles)
    {
        double_magic.push_back(-d);
    }

    for (double degrees : double_magic)
    {
        cu::mush::MoireApproximator graph_moire(tile, degrees,max_tiles);
        const auto aligned_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED);
        const auto rotated_report = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED);

        magicdeg << degrees << "    ";
        for (auto report : {aligned_report, rotated_report})
        {
            cu::mush::DeformationReport dr(report.approximation_deformation);
            magicdeg << dr.rotation_angle << "    " << dr.dilation_strain << "    " << dr.deviatoric_strain << "    ";
        }
        magicdeg << "\n";
    }
    magicdeg.close();

    return;
}

TEST_F(GrapheneTwistTest, Debug)
{
    return;
    std::ifstream magics("/home/mesto/programming/multishifter.tmp/sandbox/magic.txt");
    std::vector<double> angles;
    while(!magics.eof())
    {
        double angle;
        magics>>angle;
        angles.push_back(angle);
    }

    spit_determinants(graphene_ptr->lattice(), magic_angles, 1000);
    return

    spit_internal_error_data(graphene_ptr->lattice(), magic_angles, 3000);

    spit_internal_error_data(graphene_ptr->lattice(), magic_angles, 0);
    spit_internal_error_data(graphene_ptr->lattice(), magic_angles, 200);
    spit_internal_error_data(graphene_ptr->lattice(), magic_angles, 1000);

    spit_rotation_error_data(*graphene_ptr, magic_angles, 0);
    spit_rotation_error_data(*graphene_ptr, magic_angles, 200);
    spit_rotation_error_data(*graphene_ptr, magic_angles, 1000);
    return;

    for (int i = 0; i <= 1000; i += 100)
    {
        spit_rotation_error_data(*graphene_ptr, magic_angles, i);
    }

    double angle = 15.178178937949879; // This angle gives a coincident sqrt(3) sqrt(3) moire superlattice
    cu::mush::MoireStructureApproximator mini_graph_moire(*graphene_ptr, angle,100);
    cu::mush::MoireStructureApproximator graph_moire(*graphene_ptr, angle,100);

    const auto mini_struc_top = mini_graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED).approximate_moire_structure;
    const auto mini_struc_bottom = mini_graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED).approximate_moire_structure;
    const auto mini_stack = cu::frankenstein::stack({mini_struc_bottom, mini_struc_top});

    const auto struc_top = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED).approximate_moire_structure;
    const auto struc_bottom = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED).approximate_moire_structure;
    const auto stack = cu::frankenstein::stack({struc_bottom, struc_top});

    cu::xtal::write_poscar(mini_stack, "mini_stack.vasp");
    cu::xtal::write_poscar(stack, "scel_stack.vasp");

    for (double angle : magic_angles)
    {
        cu::mush::MoireStructureApproximator graph_moire(*graphene_ptr, angle, 1000);

        const auto a_top = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ROTATED).approximate_moire_structure;
        const auto a_bottom = graph_moire.best_smallest(ZONE::ALIGNED, LAT::ALIGNED).approximate_moire_structure;
        const auto a_stack = cu::frankenstein::stack({a_bottom, a_top});
        cu::xtal::write_poscar(a_stack, std::to_string(angle) + "aligned.vasp");

        const auto r_top = graph_moire.best_smallest(ZONE::ROTATED, LAT::ROTATED).approximate_moire_structure;
        const auto r_bottom = graph_moire.best_smallest(ZONE::ROTATED, LAT::ALIGNED).approximate_moire_structure;
        const auto r_stack = cu::frankenstein::stack({a_bottom, a_top});
        cu::xtal::write_poscar(r_stack, std::to_string(angle) + "rotated.vasp");
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
