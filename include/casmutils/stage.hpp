#include <casm/CASM_global_definitions.hh>
#include "casmutils/structure.hpp"

namespace Simplicity
{
    void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor);
    void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& MODE);
}
