#include <casmutils/xtal/lattice.hpp>

namespace casmutils
{
namespace xtal
{
Lattice::Lattice(const CASM::xtal::Lattice& init_lat) : casm_lattice(init_lat) {}
Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : casm_lattice(column_lat_mat) {}

LatticeEquals_f::LatticeEquals_f(const Lattice& ref_lat, double tol) : ref_lat(ref_lat), tol(tol) {}
bool LatticeEquals_f::operator()(const Lattice& other)
{
    return ref_lat.column_vector_matrix().isApprox(other.column_vector_matrix(), tol);
}
} // namespace xtal
} // namespace casmutils
