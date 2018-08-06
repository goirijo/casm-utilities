#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"

namespace Simplicity
{
    Rewrap::Structure make_super_structure(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat)
    {
        auto lattice_mat = struc.lattice().lat_column_mat();
        // had to cast the transformation matrix to double as Eigen does not allow mixing matrix types
        CASM::Lattice suplat(lattice_mat * col_transf_mat.cast<double>());
        return struc.create_superstruc(suplat);
    }
}
