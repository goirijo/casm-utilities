#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/strain/StrainConverter.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"

namespace Simplicity
{
    void apply_deformation_tensor(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor) 
    {
        CASM::Lattice strain_lat(deformation_tensor * struc_ptr->lattice().lat_column_mat());
        struc_ptr->set_lattice(strain_lat, CASM::FRAC);
        return;
    }

    // void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& strain_tensor, const std::string& MODE) 
    // {
    //     CASM::StrainConverter converter(MODE);
    //     auto deformation_tensor = converter.strain_metric_to_F(strain_tensor);
    //     apply_deformation_tensor(struc_ptr,deformation_tensor);
    //     return;
    // }
}
