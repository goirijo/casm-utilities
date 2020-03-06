#ifndef UTILS_LATTICE_HH
#define UTILS_LATTICE_HH

#include <casm/crystallography/Lattice.hh>
#include <casmutils/definitions.hpp>

namespace rewrap
{
class Lattice
{
public:
    Lattice() = delete;
    Lattice(const CASM::xtal::Lattice& init_lat);
    Lattice(const Eigen::Matrix3d& column_lat_mat);

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

    /// Return *this as a CASM::Lattice
    const CASM::xtal::Lattice& __get() const { return this->casm_lattice; }

private:
    CASM::xtal::Lattice casm_lattice;
};
} // namespace rewrap

namespace casmutils
{
namespace xtal
{
using rewrap::Lattice;
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

} // namespace xtal
} // namespace casmutils

#endif
