#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/crystallography/Lattice.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"

namespace Simplicity
{
    Rewrap::Structure make_super_struc(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat)
    {
        auto mylat_mat = struc.lattice().lat_column_mat();
        /* CASM::Lattice suplat(transf_mat * mylat_mat); Don't do this */
        CASM::Lattice suplat(mylat_mat * col_transf_mat.cast<double>());
        return struc.create_superstruc(suplat);
    }

    Rewrap::Structure make_super_struc(const Rewrap::Structure& struc, const Rewrap::fs::path& transf_file_path)
    {
        Eigen::Matrix3i transf_mat;
        Rewrap::fs::ifstream mat_file(transf_file_path);
        mat_file >> transf_mat;
        return make_super_struc(struc, transf_mat);
    }

}
