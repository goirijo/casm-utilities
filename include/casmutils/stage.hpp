#include <casm/CASM_global_definitions.hh>
#include "casmutils/structure.hpp"

namespace Simplicity
{
    /// Takes a pointer to a structure and applies deformation to that structure.
    /// Deforms the lattice and keeps the basis constant in fractional space.
    /// deformed_lattice = deformation_tensor * undeformed_lattice.
    /// Does not remove rotation even if it exists!!!
    void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor);

    /// Takes a pointer to a structure and applies strain to that structure.
    /// Input is unrolled strain in conventional metrics as defined in the mode.
    /// Allowed modes are 'GL' [Green-Lagrange], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky] and throws an error if the mode is not in this list.
    /// Uses functions from CASM::StrainConverter class to roll up the strain and obtain a deformation tensor.
    /// Applies deformation using apply_deformation function. 
    void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode);
}
