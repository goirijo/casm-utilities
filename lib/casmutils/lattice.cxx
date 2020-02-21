#include "casmutils/lattice.hpp"

namespace Rewrap
{
    Lattice::Lattice(const CASM::xtal::Lattice& init_lat) : CASM::xtal::Lattice(init_lat){}
    Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : CASM::xtal::Lattice(column_lat_mat){}
}
