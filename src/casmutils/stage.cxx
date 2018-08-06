#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>
#include <casm/symmetry/SymGroupRepID.hh>
#include <casm/strain/StrainConverter.hh>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/exceptions.hpp"

namespace Simplicity
{
    void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor) 
    {
        CASM::Lattice strained_lattice(deformation_tensor * struc_ptr->lattice().lat_column_mat());
        struc_ptr->set_lattice(strained_lattice, CASM::FRAC);
        return;
    }

    void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode) 
    {
        std::set<std::string> allowed_strain_metrics = {"GL", "B", "H", "EA"};
        if (allowed_strain_metrics.count(mode))
        {
            CASM::StrainConverter converter(mode);
            auto strain_tensor = converter.rollup_E(unrolled_strain);
            auto deformation_tensor = converter.strain_metric_to_F(strain_tensor);
            apply_deformation(struc_ptr, deformation_tensor);
        }
        else
        {
            throw UtilExcept::UserInputMangle(
                "Unrecognized mode. Allowed strain metrics modes are GL, B, H and EA");
        }
        return;
    }
}
