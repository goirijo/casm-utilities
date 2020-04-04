#include "casmutils/sym/cartesian.hpp"
#include <casm/crystallography/LatticeMap.hh>
#include <casmutils/mapping/structure_mapping.hpp>
#include <vector>

namespace casmutils
{
namespace mapping
{
std::vector<sym::CartOp> StructureMapper_f::make_default_point_group() const
{
    if(this->settings.use_crystal_symmetry)
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
                                     const std::vector<sym::CartOp>& init_point_group,
                                     const AllowedSpeciesType& init_allowed_species)
    : reference_structure(reference),
      lattice_to_impose(reference_structure.lattice()),
      settings(input),
      point_group(init_point_group.empty() ? make_default_point_group() : init_point_group),
      allowed_species(init_allowed_species.empty() ? make_default_allowed_species() : init_allowed_species),
      mapper(
          CASM::xtal::SimpleStrucMapCalculator(
              reference_structure.__get<CASM::xtal::SimpleStructure>(), point_group, SpecMode::ATOM, allowed_species),
          settings.strain_weight,
          settings.max_volume_change,
          settings.options,
          settings.tol,
          settings.min_va_frac,
          settings.max_va_frac)
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
    input.use_crystal_symmetry=true;
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

} // namespace mapping
} // namespace casmutils
