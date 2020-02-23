#include <casmutils/xtal/lattice.hpp>

namespace Rewrap
{
Lattice::Lattice(const CASM::xtal::Lattice& init_lat) : casm_lattice(init_lat) {}
Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : casm_lattice(column_lat_mat) {}
} // namespace Rewrap
