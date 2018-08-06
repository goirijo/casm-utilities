#include "casmutils/structure.hpp"
#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
/* #include <casm/symmetry/SymGroupRepID.hh> */
#include <casm/strain/StrainConverter.hh>
#include "casmutils/exceptions.hpp"
#include <fstream>

namespace Rewrap
{
Structure::Structure(CASM::Structure init_struc) : CASM::Structure(init_struc) {}
Structure::Structure(Rewrap::fs::path& filename) : CASM::Structure(filename) {}

Structure Structure::from_poscar(const fs::path& poscar_path)
{
    return Rewrap::Structure(CASM::Structure(poscar_path));
}

bool Structure::is_primitive() const { return CASM::Structure::is_primitive(); }

Structure Structure::primitive() const { return Simplicity::make_primitive(*this); }
} // namespace Rewrap

namespace Simplicity
{
Rewrap::Structure make_primitive(const Rewrap::Structure& input)
{
    const CASM::Structure& casted_input(input);
    CASM::Structure true_prim;
    bool is_prim = casted_input.is_primitive(true_prim);
    return true_prim;
}

Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli)
{
    CASM::Structure niggli = non_niggli;
    CASM::Lattice lat_niggli = CASM::niggli(non_niggli.lattice(), CASM::TOL);
    niggli.set_lattice(lat_niggli, CASM::CART);
    return niggli;
}

void make_niggli(Rewrap::Structure* non_niggli)
{
    CASM::Lattice lat_niggli = CASM::niggli(non_niggli->lattice(), CASM::TOL);
    non_niggli->set_lattice(lat_niggli, CASM::CART);
    return;
}

void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(printable);
    p.sort();
    p.print(outstream);
    return;
}

void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}

Rewrap::Structure make_super_structure(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat)
{
    auto lattice_mat = struc.lattice().lat_column_mat();
    // had to cast the transformation matrix to double as Eigen does not allow mixing matrix types
    CASM::Lattice suplat(lattice_mat * col_transf_mat.cast<double>());
    return struc.create_superstruc(suplat);
}

void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    CASM::Lattice strained_lattice(deformation_tensor * struc_ptr->lattice().lat_column_mat());
    struc_ptr->set_lattice(strained_lattice, CASM::FRAC);
    return;
}

Rewrap::Structure apply_deformation(const Rewrap::Structure& struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    Rewrap::Structure copy_struc(struc_ptr);
    apply_deformation(&copy_struc, deformation_tensor);
    return copy_struc;
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
        throw UtilExcept::UserInputMangle("Unrecognized mode. Allowed strain metrics modes are GL, B, H and EA");
    }
    return;
}

Rewrap::Structure apply_strain(const Rewrap::Structure& struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
{
    Rewrap::Structure copy_struc(struc_ptr);
    apply_strain(&copy_struc, unrolled_strain, mode);
    return copy_struc;
}
} // namespace Simplicity
