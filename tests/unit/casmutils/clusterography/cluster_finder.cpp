// These are classes that structure depends on
#include "../../../autotools.hh"
#include <casmutils/clusterography/cluster_finder.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>
#include <gtest/gtest.h>

namespace cu = casmutils;

class ClusterographyTest : public testing::Test
{
protected:
    void SetUp() override {}

    // Use unique pointers because Structure has no default constructor
    std::unique_ptr<casmutils::xtal::Structure> cubic_Ni_struc_ptr;
    typedef casmutils::xtal::Site Site;
    typedef casmutils::clusterography::Cluster Cluster;
};

TEST_F(ClusterographyTest, MakePeriodicOrbits)
{
    //  checks to see if Structure can be constructed
    // and has the lattice that was given associated with it.
    casmutils::fs::path pos_path(casmutils::autotools::input_filesdir / "conventional_fcc_Ni.vasp");
    casmutils::xtal::Structure pos = casmutils::xtal::Structure::from_poscar(pos_path);

    std::vector<double> max_length{3, 3, 4};
    casmutils::clusterography::make_periodic_orbits(max_length, pos);

    Eigen::Vector3d p(3, 3, 3);
    Eigen::Vector3d q(2, 2, 2);

    Site s1(p, "C");
    Site s2(q, "P");

    std::vector<Site> pq{s1, s2};
    Cluster pp(pq);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
