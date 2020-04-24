#include <algorithm>
#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Strain.hh>
#include <casm/crystallography/SuperlatticeEnumerator.hh>
#include <casm/crystallography/SymTools.hh>
#include <casm/crystallography/io/VaspIO.hh>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <fstream>
namespace
{
// surface area of a given lattice
double lattice_surface_area(const casmutils::xtal::Lattice& lat)
{
    Eigen::Vector3d a = lat[0];
    Eigen::Vector3d b = lat[1];
    Eigen::Vector3d c = lat[2];

    double ab = a.cross(b).norm();
    double bc = b.cross(c).norm();
    double ca = c.cross(b).norm();

    return std::abs(ab) + std::abs(bc) + std::abs(ca);
}

// score for determining level of boxiness, borrowed from John Goiri
double boxy_score(const casmutils::xtal::Lattice& lat)
{
    // Less surface area per volume means more boxy
    // i.e. more volume per surface area means more boxy
    return std::abs(lat.volume()) / lattice_surface_area(lat);
}
} // namespace

namespace casmutils
{
namespace xtal
{
Structure make_primitive(const Structure& input)
{
    const auto& casted_input = input.__get<CASM::xtal::BasicStructure>();
    Structure true_prim(CASM::xtal::make_primitive(casted_input));
    return true_prim;
}

Structure make_niggli(const Structure& non_niggli)
{
    Structure niggli = non_niggli;
    make_niggli(&niggli);
    return niggli;
}

void make_niggli(Structure* non_niggli)
{
    Lattice lat_niggli = make_niggli(non_niggli->lattice());
    non_niggli->set_lattice(lat_niggli, CART);
    non_niggli->within();
    return;
}

void print_poscar(const Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(printable.__get<CASM::xtal::SimpleStructure>());
    p.sort();
    p.print(outstream);
    return;
}

void write_poscar(const Structure& printable, const fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}

Structure make_superstructure(const Structure& struc, const Eigen::Matrix3i& col_transf_mat)
{
    CASM::xtal::BasicStructure superstructure(
        CASM::xtal::make_superstructure(struc.__get<CASM::xtal::BasicStructure>(), col_transf_mat));
    return Structure(superstructure);
}

void apply_deformation(Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    Lattice strained_lattice(deformation_tensor * struc_ptr->lattice().column_vector_matrix());
    struc_ptr->set_lattice(strained_lattice, FRAC);
    return;
}

Structure apply_deformation(const Structure& struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    Structure copy_struc(struc_ptr);
    apply_deformation(&copy_struc, deformation_tensor);
    return copy_struc;
}

Eigen::Matrix3d rollup_strain_metric(const Eigen::Ref<const Eigen::VectorXd>& unrolled_strain)
{
    Eigen::Matrix3d rolled_up_strain = Eigen::Matrix3d::Zero();
    rolled_up_strain(0, 0) = unrolled_strain(0);
    rolled_up_strain(1, 1) = unrolled_strain(1);
    rolled_up_strain(2, 2) = unrolled_strain(2);
    rolled_up_strain(1, 2) = unrolled_strain(3) / sqrt(2);
    rolled_up_strain(0, 2) = unrolled_strain(4) / sqrt(2);
    rolled_up_strain(0, 1) = unrolled_strain(5) / sqrt(2);
    rolled_up_strain(2, 1) = unrolled_strain(3) / sqrt(2);
    rolled_up_strain(2, 0) = unrolled_strain(4) / sqrt(2);
    rolled_up_strain(1, 0) = unrolled_strain(5) / sqrt(2);
    return rolled_up_strain;
}

void apply_strain(Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
{
    std::set<std::string> allowed_strain_metrics = {"GL", "B", "H", "EA"};
    if (allowed_strain_metrics.count(mode))
    {
        auto strain_tensor = rollup_strain_metric(unrolled_strain);
        auto deformation_tensor =
            CASM::strain::metric_to_deformation_tensor<CASM::strain::METRIC::GREEN_LAGRANGE>(strain_tensor);
        apply_deformation(struc_ptr, deformation_tensor);
    }
    else
    {
        throw except::UserInputMangle("Unrecognized mode. Allowed strain metrics modes are GL, B, H and EA");
    }
    return;
}

Structure apply_strain(const Structure& struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
{
    Structure copy_struc(struc_ptr);
    apply_strain(&copy_struc, unrolled_strain, mode);
    return copy_struc;
}

std::vector<Structure> make_superstructures_of_volume(const Structure& structure, const int volume)
{
    std::vector<Structure> all_superstructures;
    CASM::xtal::ScelEnumProps enum_props(volume, volume + 1);
    // TODO: Bring crystal group operations into casmutils
    std::vector<CASM::xtal::SymOp> pg = CASM::xtal::make_point_group(structure.lattice().__get());
    CASM::xtal::SuperlatticeEnumerator lat_enumerator(structure.lattice().__get(), pg, enum_props);

    for (const auto& lat : lat_enumerator)
    {
        Eigen::Matrix3i transfmat =
            (structure.lattice().column_vector_matrix().inverse() * lat.lat_column_mat()).cast<int>();
        Structure super = CASM::xtal::make_superstructure(structure.__get<CASM::xtal::BasicStructure>(), transfmat);
        make_niggli(&super);

        all_superstructures.emplace_back(std::move(super));
    }

    return all_superstructures;
}
} // namespace xtal
} // namespace casmutils
