// These are classes that structure depends on
#include "../../../autotools.hh"
#include <casmutils/clusterography/cluster_finder.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <functional>
#include <gtest/gtest.h>

namespace cu = casmutils;

// useful function for site filter
bool site_filter(const cu::xtal::Site& site) { return true; }

class ClusterographyTest : public testing::Test
{
protected:
    void SetUp() override
    {
        cu::fs::path pos_path(cu::autotools::input_filesdir / "conventional_fcc_Ni.vasp");
        struc_ptr.reset(new cu::xtal::Structure(cu::xtal::Structure::from_poscar(pos_path)));
        max_length0 = {3, 2, 4};
        max_length1 = {3, 2, 4.01};
        cu_site_filter = site_filter;
        tol = 1e-5;
    }

    std::function<bool(cu::xtal::Site)> cu_site_filter;
    std::unique_ptr<casmutils::xtal::Structure> struc_ptr;
    std::vector<double> max_length0;
    std::vector<double> max_length1;
    double tol;
    typedef casmutils::xtal::Site Site;
    typedef casmutils::clusterography::Cluster Cluster;
    typedef casmutils::clusterography::Orbit Orbit;
};

TEST_F(ClusterographyTest, MakePeriodicOrbits)
{
    std::vector<Orbit> orbits0 = cu::clusterography::make_periodic_orbits(max_length0, *struc_ptr);
    std::vector<Orbit> orbits1 = cu::clusterography::make_periodic_orbits(max_length1, *struc_ptr);

    // In the first case, you will have 3 orbits (1 null point, 1 single point, 1 two point clusters)
    // There will be 4 clusters in single point cluster orbit (4 atoms in the conventional fcc)
    // There will be 24 clusters in two point cluster orbit (face centered & corner atom cluster - 4 clusters
    // for face - 6 faces - 24 clusters in orbit)
    // This illustrates that function returns clusters which are less than max_length for a branch (not equal
    // to max_length)
    EXPECT_EQ(orbits0.size(), 3);
    EXPECT_EQ(orbits0[1].size(), 4);
    EXPECT_EQ(orbits0[2].size(), 24);

    // In the second case, you will have 4 orbits (1 null point, 1 single point, 2 two point clusters)
    // The second 2 point cluster orbit now includes corner - corner atoms in the unit cell
    EXPECT_EQ(orbits1.size(), 4);
    EXPECT_EQ(orbits1[1].size(), 4);
    EXPECT_EQ(orbits1[2].size(), 24);
    EXPECT_EQ(orbits1[3].size(), 12);
}

TEST_F(ClusterographyTest, MakeLocalOrbits)
{
    std::vector<double> cut_off_radius{1, 3, 3};

    // Phenomenal Cluster - just a single point cluster with the site at origin
    Cluster phenomenal_cluster;
    phenomenal_cluster.add_site(struc_ptr->basis_sites()[0]);

    std::vector<Orbit> orbits0 =
        cu::clusterography::make_local_orbits(phenomenal_cluster, false, max_length0, cut_off_radius, *struc_ptr);

    // When you don't include the phenomenal cluster, there will be 12 sites around the origin
    EXPECT_EQ(orbits0.size(), 3);
    EXPECT_EQ(orbits0[1].size(), 12);
    EXPECT_EQ(orbits0[2].size(), 24);

    std::vector<Orbit> orbits1 =
        cu::clusterography::make_local_orbits(phenomenal_cluster, true, max_length0, cut_off_radius, *struc_ptr);

    // When you include the phenomenal cluster 2 new orbits will be generated
    // orbit of branch 1 with origin as a cluster
    // orbit of branch 2 with origin & face centered atoms as a cluster (which will be 12)
    EXPECT_EQ(orbits1.size(), 5);
    EXPECT_EQ(orbits1[1].size(), 1);
    EXPECT_EQ(orbits1[2].size(), 12);
    EXPECT_EQ(orbits1[3].size(), 12);
    EXPECT_EQ(orbits1[4].size(), 24);
}

TEST_F(ClusterographyTest, ApplySymOpToCluster)
{
    Cluster test_cluster;
    test_cluster.add_site(struc_ptr->basis_sites()[0]);
    test_cluster.add_site(struc_ptr->basis_sites()[1]);

    Eigen::Matrix3d r90;
    r90 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    cu::sym::CartOp rotation_90 = cu::sym::CartOp::point_operation(r90);

    Cluster new_cluster = rotation_90 * test_cluster;
    EXPECT_EQ(new_cluster.size(), 2);
    EXPECT_TRUE(new_cluster[0].cart().isZero(tol));

    Eigen::Vector3d new_site(-2, 2, 0);
    EXPECT_TRUE(new_cluster[1].cart().isApprox(new_site));
}

// TODO: Add more tests for cluster & orbit functions
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
