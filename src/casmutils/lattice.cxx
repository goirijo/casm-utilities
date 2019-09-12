#include "casmutils/lattice.hpp"

namespace Rewrap
{
    Lattice::Lattice(const CASM::Lattice& init_lat) : CASM::Lattice(init_lat){}
    Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : CASM::Lattice(column_lat_mat){}
}
