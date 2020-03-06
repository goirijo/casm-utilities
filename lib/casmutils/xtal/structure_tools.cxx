#include <algorithm>
#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/LatticeMap.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/SuperlatticeEnumerator.hh>
#include <casm/crystallography/SymTools.hh>
#include <casm/crystallography/io/VaspIO.hh>
#include <casm/strain/StrainConverter.hh>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/stage.hpp>
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
    Lattice lat_niggli =
        CASM::xtal::niggli(CASM::xtal::Lattice(non_niggli->lattice().column_vector_matrix()), CASM::TOL);
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

Structure make_super_structure(const Structure& struc, const Eigen::Matrix3i& col_transf_mat)
{
    auto lattice_mat = struc.lattice().column_vector_matrix();
    // had to cast the transformation matrix to double as Eigen does not allow mixing matrix types
    CASM::xtal::Lattice suplat(lattice_mat * col_transf_mat.cast<double>());
    return Structure(struc.__get<CASM::xtal::BasicStructure>().create_superstruc(suplat));
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

void apply_strain(Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
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

mapping::MappingNode structure_map(const Structure& map_reference_struc, const Structure& mappable_struc)
{
    throw except::NotImplemented();
}

std::pair<double, double> structure_score(const mapping::MappingNode& mapping_data)
{
    double lattice_score = CASM::xtal::StrainCostCalculator::iso_strain_cost(
        mapping_data.stretch,
        mapping_data.child.column_vector_matrix().determinant() / std::max(int(mapping_data.permutation.size()), 1));
    double basis_score = (mapping_data.stretch.inverse() * mapping_data.displacement).squaredNorm() /
                         double(std::max(int(mapping_data.permutation.size()), 1));
    return std::make_pair(lattice_score, basis_score);
}

mapping::MappingNode structure_map(const mapping::MappingInput& input, const Structure& mappable_struc)
{
    std::vector<std::unordered_set<std::string>> allowed_at_sites(input.parent.basis_sites().size(),
                                                                  input.allowed_species);
    CASM::xtal::SimpleStrucMapCalculator calc_interface(input.parent.__get<CASM::xtal::SimpleStructure>(),
                                                        input.point_group, input.mode, allowed_at_sites);
    CASM::xtal::StrucMapper mapper(calc_interface, input.strain_weight, input.max_volume_change, input.options,
                                   input.tol, input.min_va_frac, input.max_va_frac);
    if (input.is_ideal)
    {
        return mapper.map_ideal_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), input.num_best_maps);
    }

    return *(mapper
                 .map_deformed_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), input.num_best_maps,
                                     input.max_cost, input.min_cost, input.keep_invalid_mapping_nodes)
                 .begin());
}

std::vector<Structure> make_superstructures_of_volume(const Structure& structure, const int volume)
{
    std::vector<Structure> all_superstructures;
    CASM::xtal::ScelEnumProps enum_props(volume, volume + 1);
    std::vector<CASM::xtal::SymOp> pg = CASM::xtal::make_point_group(structure.lattice().__get());
    CASM::xtal::SuperlatticeEnumerator lat_enumerator(structure.lattice().__get(), pg, enum_props);

    for (const auto& lat : lat_enumerator)
    {
        Structure super = structure.__get<CASM::xtal::BasicStructure>().create_superstruc(lat);
        make_niggli(&super);

        all_superstructures.emplace_back(std::move(super));
    }

    return all_superstructures;
}
} // namespace xtal
} // namespace casmutils
