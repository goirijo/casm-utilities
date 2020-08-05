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

/// This functor class provides a unary predicate equals function
/// for casmutils::xtal::Lattice .
class LatticeEquals_f
{
public:
    /// The comparator requires a tolerance
    LatticeEquals_f(const Lattice& ref_lat, double tol);
    /// returns true is ref_lat is equal to other
    bool operator()(const Lattice& other);

private:
    Lattice ref_lat;
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
