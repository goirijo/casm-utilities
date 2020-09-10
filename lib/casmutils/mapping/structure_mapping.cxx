#include "casmutils/sym/cartesian.hpp"
#include <casm/crystallography/LatticeMap.hh>
#include <casmutils/mapping/structure_mapping.hpp>
#include <vector>

#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/symmetry.hpp>
namespace casmutils
{
namespace mapping
{

namespace
{
double atomic_cost_child(const MappingReport& mapped_result, int Nsites)
{
    Nsites = std::max(Nsites, int(1));
    double atomic_vol = mapped_result.reference_lattice.volume() / double(Nsites) / mapped_result.stretch.determinant();
    return pow(3. * abs(atomic_vol) / (4. * M_PI), -2. / 3.) *
           (mapped_result.stretch.inverse() * mapped_result.displacement).squaredNorm() / double(Nsites);
}
//*******************************************************************************************

double atomic_cost_parent(const MappingReport& mapped_result, int Nsites)
{
    Nsites = std::max(Nsites, int(1));
    // mean square displacement distance in deformed coordinate system
    double atomic_vol = mapped_result.reference_lattice.volume() / double(Nsites);
    return pow(3. * abs(atomic_vol) / (4. * M_PI), -2. / 3.) * (mapped_result.displacement).squaredNorm() /
           double(Nsites);
}

//*******************************************************************************************

double atomic_cost(const MappingReport& mapped_result, int Nsites)
{
    // mean square displacement distance in deformed coordinate system
    return (atomic_cost_child(mapped_result, Nsites) + atomic_cost_parent(mapped_result, Nsites)) / 2.;
}

} // namespace

std::vector<sym::CartOp> StructureMapper_f::make_default_factor_group() const
{
    if (this->settings.use_crystal_symmetry)
    {
        return xtal::make_factor_group(reference_structure,this->settings.tol);
    }

    return {sym::CartOp::identity()};
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
    double lattice_score = CASM::xtal::StrainCostCalculator::isotropic_strain_cost(mapping_data.stretch);
    double basis_score = atomic_cost(mapping_data, std::max(int(mapping_data.permutation.size()), 1));
    return std::make_pair(lattice_score, basis_score);
}

} // namespace mapping
} // namespace casmutils
