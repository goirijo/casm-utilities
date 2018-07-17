#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include "casmutils/structure.hpp"

namespace Simplicity
{
    Rewrap::Structure make_super_struc(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat);
    Rewrap::Structure make_super_struc(const Rewrap::Structure& struc, const Rewrap::fs::path& transf_file_path);
}
