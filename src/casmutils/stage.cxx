#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/symmetry/SymGroupRepID.hh>
#include <casm/strain/StrainConverter.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"

namespace Simplicity
{
    // Applies deformation to a structure
    // deforms the lattice and keeps the basis constant in fractional space
    // deformed_lattice = deformation_tensor * undeformed_lattice
    // even if exists does not remove rotation!!!
    void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor) 
    {
        CASM::Lattice strain_lat(deformation_tensor * struc_ptr->lattice().lat_column_mat());
        struc_ptr->set_lattice(strain_lat, CASM::FRAC);
        return;
    }

    // Applies strain to a structure. Input is unrolled strain in conventional metrics as defined in the MODE
    // Uses functions from CASM::StrainConverter class to roll up the strain and obtain a deformation tensor
    // applies deformation using apply_deformation function. 
    void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& MODE) 
    {
        CASM::StrainConverter converter(MODE);
        auto strain_tensor = converter.rollup_E(unrolled_strain);
        auto deformation_tensor = converter.strain_metric_to_F(strain_tensor);
        apply_deformation(struc_ptr, deformation_tensor);
        return;
    }
}
