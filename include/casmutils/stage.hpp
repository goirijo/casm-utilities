#include <casm/CASM_global_definitions.hh>
#include "casmutils/structure.hpp"

namespace Simplicity
{
    /// Returns a super structure after applying a transformation matrix to the structure.
    /// Applies transformation to the lattice and uses CASM::Structure::create_superstruc to fill the basis. 
    /// transformed_lattice =  original_lattice * transformation_matrix
    Rewrap::Structure make_super_structure(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat);
}
