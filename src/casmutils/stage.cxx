#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"

namespace Simplicity
{
    // returns a structure after applying transformed matrix to the structure
    // apples transformation to the lattice and uses CASM::Structure::create_superstruc to fill the basis 
    // transformed_lattice =  original_lattice * transformation_matrix
    Rewrap::Structure make_super_structure(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat)
    {
        auto mylat_mat = struc.lattice().lat_column_mat();
        // had to cast the transformation matrix to double as Eigen does not allow mixing matrix types
        CASM::Lattice suplat(mylat_mat * col_transf_mat.cast<double>());
        return struc.create_superstruc(suplat);
    }

}
