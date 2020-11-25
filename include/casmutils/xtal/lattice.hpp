#ifndef UTILS_LATTICE_HH
#define UTILS_LATTICE_HH

#include <casm/crystallography/Lattice.hh>
#include <casmutils/definitions.hpp>

namespace casmutils
{
namespace xtal
{
class Lattice
{
public:
    Lattice() = delete;
    Lattice(const CASM::xtal::Lattice& init_lat);
    Lattice(const Eigen::Matrix3d& column_lat_mat);

    /// Construct lattice by specifying each individual vector
    Lattice(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

    static Lattice from_column_vector_matrix(const Eigen::Matrix3d& column_vector_matrix)
    {
        return Lattice(column_vector_matrix);
    }

    static Lattice from_row_vector_matrix(const Eigen::Matrix3d& row_vector_matrix)
    {
        return Lattice(row_vector_matrix.transpose());
    }

    /// Constructing lattice by specifying lattice parameters
    /// Since there can be many lattices which can be constructed using the given lattice paramters,
    /// the following asssumptions are made:
    /// z-axis is aligned with "c" vector and "a" vector lies in the x-z plane
    /// The values of angles should be provided in degrees
    /// alpha - angle between b & c
    /// beta - angle between c & a
    /// gamma - angle between a & b
    static Lattice from_lattice_parameters(const double& a,
                                           const double& b,
                                           const double& c,
                                           const double& alpha,
                                           const double& beta,
                                           const double& gamma);

    // TODO: Read up on Eigen::Matrix3d::ColXpr and decide if you prefer this. CASM does it this way.
    /// Return the ith vector of the lattice
    Eigen::Vector3d operator[](int i) const { return this->__get()[i]; }

    /// Return the first vector of the lattice
    Eigen::Vector3d a() const { return this->operator[](0); }

    /// Return the second vector of the lattice
    Eigen::Vector3d b() const { return this->operator[](1); }

    /// Return the third vector of the lattice
    Eigen::Vector3d c() const { return this->operator[](2); }

    /// Return the volume of this lattice
    double volume() const { return this->__get().volume(); }

    /// Returns the matrix representation of this lattice where each lattice vector is a column
    Eigen::Matrix3d column_vector_matrix() const { return this->__get().lat_column_mat(); }

    /// Returns the matrix representation of this lattice where each lattice vector is a row
    Eigen::Matrix3d row_vector_matrix() const { return this->column_vector_matrix().transpose(); }

    /// Return *this as a CASM::Lattice
    const CASM::xtal::Lattice& __get() const { return this->casm_lattice; }

private:
    static Eigen::Matrix3d
    stack_column_vectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c);

    CASM::xtal::Lattice casm_lattice;
};

// TODO: Make this a binary comparator, fix is_equal, and implement UnaryComparator_f
/// This functor class provides a unary predicate equals function
/// for casmutils::xtal::Lattice .
class LatticeEquals_f
{
public:
    /// The comparator requires a tolerance
    LatticeEquals_f(double tol);
    /// returns true is ref_lat is equal to other by direct vector comparison
    bool operator()(const Lattice& ref_lat, const Lattice& other) const;

private:
    double tol;
};

/// True if the two lattices are related by a unimodular transformation
class LatticeIsEquivalent_f
{
public:
    LatticeIsEquivalent_f(double tol) : tol(tol) {}
    /// True if the two lattices are related by a unimodular transformation
    /// (equivalent under point group operation)
    bool operator()(const Lattice& reference, const Lattice& other) const;

private:
    double tol;
};

void make_niggli(Lattice* lattice_ptr);
Lattice make_niggli(const Lattice& non_niggli_lattice);

/// Retrun the reciprocal of the given lattice, where each reciprocal vector
/// is perpendicular to the other two real ones
Lattice make_reciprocal(const Lattice& real_lattice);

/// For a deformation matrix F, return the polar decomposition of its rotational and
/// strain components R (first) and U (second), where F=R*U
std::pair<Eigen::Matrix3d, Eigen::Matrix3d> polar_decomposition(Eigen::Matrix3d const& F);

/// Invert the directions of the lattice vectors if necessary, such that
/// the column vector matrix has a positive determinant
Lattice make_right_handed(const Lattice& left_handed_lattice);

/// Create a superlattice using the provided integer transformation matrix
xtal::Lattice make_superlattice(const xtal::Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat);

/// Given a lattice and a vector of integer Miller indices, return the smallest superlattice
/// that has the a and b vectors spanning the specified plane
xtal::Lattice slice_along_plane(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);
} // namespace xtal
} // namespace casmutils

#endif
