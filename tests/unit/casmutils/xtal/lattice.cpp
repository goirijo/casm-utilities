#include <casmutils/misc.hpp>
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

        fcc_ptr.reset(new Lattice(fcc_matrix));
        fcc_copy_ptr.reset(new Lattice(fcc_matrix));

        Eigen::Matrix3d conventional_fcc_matrix;
        conventional_fcc_matrix << 3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0;

        conventional_fcc_ptr.reset(new Lattice(conventional_fcc_matrix));

        Eigen::Matrix3d non_niggli_conventional_fcc_matrix;
        non_niggli_conventional_fcc_matrix << 3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 3.0, 3.0, 3.0;

        non_niggli_conventional_fcc_ptr.reset(new Lattice(non_niggli_conventional_fcc_matrix));
        Eigen::Matrix3d bcc_matrix;
        bcc_matrix << -1.5, 1.5, 1.5, 1.5, -1.5, 1.5, 1.5, 1.5, -1.5;
        bcc_ptr.reset(new Lattice(bcc_matrix));

        Eigen::Matrix3d hcp_matrix;
        hcp_matrix << -2.871, -1.4355, 0.0, 0.0, 2.486358934265, 0.0, 0.0, 0.0, 4.635;
        hcp_ptr.reset(new Lattice(hcp_matrix));

        lat_by_vector.reset(new Lattice(a, b, c));
    }
    using Lattice = casmutils::xtal::Lattice;
    // Use unique pointers because Lattice has no default constructor
    std::unique_ptr<Lattice> fcc_ptr;
    std::unique_ptr<Lattice> fcc_copy_ptr;

    std::unique_ptr<Lattice> conventional_fcc_ptr;
    std::unique_ptr<Lattice> non_niggli_conventional_fcc_ptr;

    std::unique_ptr<Lattice> bcc_ptr;
    std::unique_ptr<Lattice> hcp_ptr;

    Eigen::Vector3d a = Eigen::Vector3d(1, 1, 1);
    Eigen::Vector3d b = Eigen::Vector3d(2, 2, 2);
    Eigen::Vector3d c = Eigen::Vector3d(3, 3, 3);
    std::unique_ptr<Lattice> lat_by_vector;
    double tol = 1e-5;
};

TEST_F(LatticeTest, ContructByVector)
{
    // Construct the lattice one vector at a time
    EXPECT_TRUE(lat_by_vector->a() == a);
    EXPECT_TRUE(lat_by_vector->b() == b);
    EXPECT_TRUE(lat_by_vector->c() == c);
}

TEST_F(LatticeTest, ColumnMatrixCheck)
{
    Eigen::Matrix3d column_matrix = lat_by_vector->column_vector_matrix();
    EXPECT_TRUE(Eigen::Vector3d(column_matrix.col(0)) == a);
    EXPECT_TRUE(Eigen::Vector3d(column_matrix.col(1)) == b);
    EXPECT_TRUE(Eigen::Vector3d(column_matrix.col(2)) == c);
}

TEST_F(LatticeTest, RowMatrixCheck)
{
    Eigen::Matrix3d row_matrix = lat_by_vector->row_vector_matrix();
    EXPECT_TRUE(Eigen::Vector3d(row_matrix.row(0)) == a);
    EXPECT_TRUE(Eigen::Vector3d(row_matrix.row(1)) == b);
    EXPECT_TRUE(Eigen::Vector3d(row_matrix.row(2)) == c);
}

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
    Eigen::Vector3d unitvec = (*hcp_ptr)[2] / (*hcp_ptr)[2].norm();
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
    casmutils::xtal::LatticeEquals_f is_equal_to_fcc_lattice(*fcc_ptr, tol);
    EXPECT_TRUE(is_equal_to_fcc_lattice(*fcc_copy_ptr));
    EXPECT_FALSE(is_equal_to_fcc_lattice(fcc_with_distortion_lat));
}

TEST_F(LatticeTest, MakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    casmutils::xtal::make_niggli(non_niggli_conventional_fcc_ptr.get());
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(
        *conventional_fcc_ptr, *non_niggli_conventional_fcc_ptr, tol));
}

TEST_F(LatticeTest, ConstMakeNiggli)
{
    // checks to see if you can make a skewed cell as niggli
    const Lattice niggli = casmutils::xtal::make_niggli(*non_niggli_conventional_fcc_ptr);
    EXPECT_TRUE(casmutils::is_equal<casmutils::xtal::LatticeEquals_f>(*conventional_fcc_ptr, niggli, tol));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
