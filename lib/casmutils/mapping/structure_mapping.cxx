#include "casmutils/sym/cartesian.hpp"
#include <casm/crystallography/LatticeMap.hh>
#include <casmutils/mapping/structure_mapping.hpp>
#include <vector>

namespace casmutils
{
namespace mapping
{
std::vector<sym::CartOp> StructureMapper_f::make_default_factor_group() const
{
    if (this->settings.use_crystal_symmetry)
    {
        return CASM::xtal::make_factor_group(this->reference_structure.__get<CASM::xtal::BasicStructure>());
    }

    return std::vector<sym::CartOp>{sym::CartOp::identity()};
}

StructureMapper_f::AllowedSpeciesType StructureMapper_f::make_default_allowed_species() const
{
    AllowedSpeciesType default_allowed_species;
    for (const auto& site : reference_structure.basis_sites())
    {
        std::vector<std::string> at_site_occs;
        at_site_occs.push_back(site.label());
        default_allowed_species.push_back(at_site_occs);
    }
    return default_allowed_species;
}

typedef CASM::xtal::SimpleStructure::SpeciesMode SpecMode;
StructureMapper_f::StructureMapper_f(const xtal::Structure& reference,
                                     const MappingInput& input,
                                     const std::vector<sym::CartOp>& init_factor_group,
                                     const AllowedSpeciesType& init_allowed_species)
    : reference_structure(reference),
      lattice_to_impose(reference_structure.lattice()),
      settings(input),
      factor_group(init_factor_group.empty() ? make_default_factor_group() : init_factor_group),
      allowed_species(init_allowed_species.empty() ? make_default_allowed_species() : init_allowed_species),
      mapper(
          CASM::xtal::SimpleStrucMapCalculator(
              reference_structure.__get<CASM::xtal::SimpleStructure>(), factor_group, SpecMode::ATOM, allowed_species),
          settings.strain_weight,
          settings.max_volume_change,
          settings.options,
          settings.tol,
          settings.min_vacancy_fraction,
          settings.max_vacancy_fraction)
{
    // Apologies for the ugly constructor we need to unpack input into
    // its individual values and do some layered inline construction
    //
    // explain more pls. what is "layered inline construction"?
}

std::vector<MappingReport> StructureMapper_f::operator()(const xtal::Structure& mappable_struc) const
{
    if (settings.assume_ideal_structure)
    {
        return this->ideal_map(mappable_struc);
    }

    return this->map(mappable_struc);
}

std::vector<MappingReport> StructureMapper_f::map(const xtal::Structure& mappable_struc) const
{

    auto casmnodes = mapper.map_deformed_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(),
                                               settings.k_best_maps,
                                               settings.max_cost,
                                               settings.min_cost,
                                               settings.keep_invalid_mapping_nodes);
    std::vector<MappingReport> casted_set(casmnodes.begin(), casmnodes.end());
    return casted_set;
}

std::vector<MappingReport> StructureMapper_f::ideal_map(const xtal::Structure& mappable_struc) const
{
    auto casmnodes = mapper.map_ideal_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(),
                                            settings.k_best_maps,
                                            settings.max_cost,
                                            settings.min_cost,
                                            settings.keep_invalid_mapping_nodes);
    std::vector<MappingReport> casted_set(casmnodes.begin(), casmnodes.end());
    return casted_set;
}

//***********************************************************************************//

std::vector<mapping::MappingReport> map_structure(const xtal::Structure& map_reference_struc,
                                                  const xtal::Structure& mappable_struc)
{
    mapping::MappingInput input;
    input.use_crystal_symmetry = true;
    mapping::StructureMapper_f mapper(map_reference_struc, input);
    return mapper(mappable_struc);
}

std::pair<double, double> structure_score(const mapping::MappingReport& mapping_data)
{
    double lattice_score = CASM::xtal::StrainCostCalculator::iso_strain_cost(
        mapping_data.stretch,
        mapping_data.mapped_lattice.column_vector_matrix().determinant() /
            std::max(int(mapping_data.permutation.size()), 1));
    double basis_score = (mapping_data.stretch.inverse() * mapping_data.displacement).squaredNorm() /
                         double(std::max(int(mapping_data.permutation.size()), 1));
    return std::make_pair(lattice_score, basis_score);
}

mapping::MappingReport symmetry_preserving_mapping_report(const mapping::MappingReport& mapping_data,
                                                          const std::vector<sym::CartOp>& group_as_operations,
                                                          const std::vector<sym::PermRep>& group_as_permutations)
{

    const auto disp_matrix = mapping_data.displacement;
    auto symmetry_preserving_displacement = disp_matrix;
    auto symmetry_preserving_stretch = mapping_data.stretch;
    symmetry_preserving_displacement.setZero();
    symmetry_preserving_stretch.setZero();
    for (int i = 0; i < group_as_operations.size(); ++i)
    {
        auto transformed_disp = group_as_operations[i].matrix * disp_matrix;
        Eigen::MatrixXd transformed_and_permuted_disp = transformed_disp;
        transformed_and_permuted_disp.setZero();
        int ind = 0;
        for (const auto& j : group_as_permutations[i])
        {
            transformed_and_permuted_disp.col(j) += transformed_disp.col(ind);
            ind++;
        }
        symmetry_preserving_displacement += transformed_and_permuted_disp / group_as_operations.size();
        Eigen::MatrixXd transformed_stretch =
            group_as_operations[i].matrix.transpose() * mapping_data.stretch * group_as_operations[i].matrix;
        symmetry_preserving_stretch += transformed_stretch / group_as_operations.size();
    }
    auto new_report = mapping_data;
    new_report.stretch = symmetry_preserving_stretch;
    new_report.displacement = symmetry_preserving_displacement;
    return new_report;
}
} // namespace mapping
} // namespace casmutils
