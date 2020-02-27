#include <casmutils/xtal/lattice.hpp>
#include <gtest/gtest.h>
#include <memory>

class LatticeTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Matrix3d fcc_matrix;
        fcc_matrix << 0.0, 1.5, 1.5, 1.5, 0.0, 1.5, 1.5, 1.5, 0.0;

        fcc_ptr.reset(new rewrap::Lattice(fcc_matrix));
        fcc_copy_ptr.reset(new rewrap::Lattice(fcc_matrix));

        Eigen::Matrix3d bcc_matrix;
        bcc_matrix << -1.5, 1.5, 1.5, 1.5, -1.5, 1.5, 1.5, 1.5, -1.5;
        bcc_ptr.reset(new rewrap::Lattice(bcc_matrix));

        Eigen::Matrix3d hcp_matrix;
        hcp_matrix << -2.871, -1.4355, 0.0, 0.0, 2.486358934265, 0.0, 0.0, 0.0, 4.635;
        hcp_ptr.reset(new rewrap::Lattice(hcp_matrix));
    }

    // Use unique pointers because Lattice has no default constructor
    std::unique_ptr<rewrap::Lattice> fcc_ptr;
    std::unique_ptr<rewrap::Lattice> fcc_copy_ptr;

    std::unique_ptr<rewrap::Lattice> bcc_ptr;

    std::unique_ptr<rewrap::Lattice> hcp_ptr;
};

TEST_F(LatticeTest, ConstructandGetMatrix)
{
    //  constructs two different objects with the same matrix
    //  examines the functionality of the constructor and
    // column_vector_matrix()
    EXPECT_EQ(fcc_ptr->column_vector_matrix(), fcc_copy_ptr->column_vector_matrix());
}

TEST_F(LatticeTest, BracketOperator)
{
    //  uses the bracket operator to check that the
    // third vector of z oriented hcp is along the z direction
    auto unitvec = (*hcp_ptr)[2] / (*hcp_ptr)[2].norm();
    EXPECT_EQ(unitvec, Eigen::Vector3d(0, 0, 1));
}

TEST_F(LatticeTest, VectorAccess)
{
    //  grabs the vectors of a bcc lattice and sees that
    // the lengths are the same. Does the same for a and b of hcp
    // tests a(), b(), c()
    EXPECT_EQ(bcc_ptr->a().norm(), bcc_ptr->b().norm());
    EXPECT_EQ(bcc_ptr->a().norm(), bcc_ptr->c().norm());
    EXPECT_TRUE((hcp_ptr->a().norm() - hcp_ptr->b().norm()) < 1e-5);
}
TEST_F(LatticeTest, LatticeEquals)
{
    Eigen::Matrix3d fcc_with_distortion_matrix;
    fcc_with_distortion_matrix << 0.0, 1.5, 1.5, 1.5, 0.0, 1.5, 1.5, 1.5 + 1e-4, 0.0;
    casmutils::xtal::Lattice fcc_with_distortion_lat(fcc_with_distortion_matrix);
    //  checks the ability to determine if two lattices
    // are equivalent within numerical tolerance
    double tol = 1e-5;
    casmutils::xtal::LatticeEquals_f is_equal_to_fcc_lattice(*fcc_ptr, tol);
    EXPECT_TRUE(is_equal_to_fcc_lattice(*fcc_copy_ptr));
    EXPECT_FALSE(is_equal_to_fcc_lattice(fcc_with_distortion_lat));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
