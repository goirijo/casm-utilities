#include "casmutils/stage.hpp"
namespace casmutils
{
namespace mapping
{

StructureMapper::StructureMapper(const MappingInput& input)
    : mapper(CASM::xtal::SimpleStrucMapCalculator(input.parent.__get<CASM::xtal::SimpleStructure>(), input.point_group,
                                                  input.mode, input.allowed_species),
             input.strain_weight, input.max_volume_change, input.options, input.tol, input.min_va_frac,
             input.max_va_frac),
      settings(input)
{
    // Apologies for the ugly constructor we need to unpack input into
    // its individual values and do some layered inline construction
}

std::vector<MappingNode> StructureMapper::map(const xtal::Structure& mappable_struc) const
{

    auto casmnodes =
        mapper.map_deformed_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), settings.num_best_maps,
                                  settings.max_cost, settings.min_cost, settings.keep_invalid_mapping_nodes);
    std::vector<MappingNode> casted_set(casmnodes.begin(), casmnodes.end());
    return casted_set;
}
std::vector<MappingNode> StructureMapper::ideal_map(const xtal::Structure& mappable_struc) const
{
    auto casmnodes = mapper.map_ideal_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), settings.num_best_maps,
                                            settings.max_cost, settings.min_cost, settings.keep_invalid_mapping_nodes);
    std::vector<MappingNode> casted_set(casmnodes.begin(), casmnodes.end());
    return casted_set;
}
} // namespace mapping
} // namespace casmutils
