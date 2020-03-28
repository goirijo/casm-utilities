#include "casmutils/stage.hpp"

namespace casmutils
{
namespace mapping
{

typedef CASM::xtal::SimpleStructure::SpeciesMode SpecMode;
StructureMapper_f::StructureMapper_f(const MappingInput& input)
    : mapper(CASM::xtal::SimpleStrucMapCalculator(input.reference_structure.__get<CASM::xtal::SimpleStructure>(),
                                                  input.point_group,
                                                  SpecMode::ATOM,
                                                  input.allowed_species),
             input.strain_weight,
             input.max_volume_change,
             input.options,
             input.tol,
             input.min_va_frac,
             input.max_va_frac),
      settings(input)
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
} // namespace mapping
} // namespace casmutils
