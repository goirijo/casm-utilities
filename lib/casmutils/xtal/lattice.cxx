#include <casmutils/xtal/lattice.hpp>

namespace casmutils
{
namespace xtal
{
Lattice::Lattice(const CASM::xtal::Lattice& init_lat) : casm_lattice(init_lat) {}
Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : casm_lattice(column_lat_mat) {}

Lattice::Lattice(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
    : Lattice(Lattice::stack_column_vectors(a, b, c))
{
}

Eigen::Matrix3d Lattice::stack_column_vectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                              const Eigen::Vector3d& c)
{
    Eigen::Matrix3d column_matrix;
    column_matrix.col(0) = a;
    column_matrix.col(1) = b;
    column_matrix.col(2) = c;
    return column_matrix;
}

LatticeEquals_f::LatticeEquals_f(const Lattice& ref_lat, double tol) : ref_lat(ref_lat), tol(tol) {}
bool LatticeEquals_f::operator()(const Lattice& other)
{
    return ref_lat.column_vector_matrix().isApprox(other.column_vector_matrix(), tol);
}
} // namespace xtal
} // namespace casmutils
