#include "casmutils/stage.hpp"
namespace casmutils
{
namespace mapping
{

StructureMapper::StructureMapper(const MappingInput& input)
    : mapper(CASM::xtal::SimpleStrucMapCalculator(input.parent.__get<CASM::xtal::SimpleStructure>(), input.point_group,
                                                  input.mode, input.allowed_species),
             input.strain_weight, input.max_volume_change, input.options, input.tol, input.min_va_frac,
             input.max_va_frac)
{
    // Apologies for the ugly constructor we need to unpack input into
    // its individual values and do some layered inline construction
}

MappingNode StructureMapper::map(const xtal::Structure& mappable_struc) const
{

    // This call populates default values for number of maps to 1
    // maximum mapping cost to 1e9
    // minimum mapping cost to 0
    // and it keeps invalid mapping nodes so we can check later
    return *(mapper.map_deformed_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), 1, 1e9, 0, true).begin());
}
MappingNode StructureMapper::ideal_map(const xtal::Structure& mappable_struc) const
{

    // This call populates default values for number of maps to 1
    // maximum mapping cost to 1e9
    // minimum mapping cost to 0
    // and it keeps invalid mapping nodes so we can check later
    return *(mapper.map_ideal_struc(mappable_struc.__get<CASM::xtal::SimpleStructure>(), 1, 1e9, 0, true).begin());
}
} // namespace mapping
} // namespace casmutils
