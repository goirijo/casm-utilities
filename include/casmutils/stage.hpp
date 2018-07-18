#include <casm/CASM_global_definitions.hh>
#include "casmutils/structure.hpp"

namespace Simplicity
{
    void apply_deformation_tensor(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor);
    void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& strain_tensor, const std::string& MODE);
}
