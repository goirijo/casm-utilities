#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casmutils/exceptions.hpp"
#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/structure.hpp"
namespace casmutils
{
namespace mapping
{

/// Holds the results of a structure map
/// Fundamentally this includes some strain representation,
/// displacement representation, site assignment matrix/permutation vector,
/// parent superlattice, rigid rotation, and rigid translation.
/// These quantities can be transformed and then converted to a mapping score
/// using strain_cost and basis_cost
struct MappingNode
{
    MappingNode(const CASM::xtal::MappingNode& casm_mapping_node)
        : isometry(casm_mapping_node.isometry()),
          stretch(casm_mapping_node.stretch()),
          translation(casm_mapping_node.translation()),
          displacement(casm_mapping_node.atom_displacement),
          parent(casm_mapping_node.lat_node.parent.superlattice()),
          child(casm_mapping_node.lat_node.child.superlattice())
    {
        std::vector<int> p(casm_mapping_node.atom_permutation.begin(), casm_mapping_node.atom_permutation.end());
        permutation = p;
        if (!casm_mapping_node.is_valid)
        {
            throw except::InvalidMap();
        }
    }
    Eigen::Matrix3d isometry;
    Eigen::Matrix3d stretch;
    Eigen::Vector3d translation;
    Eigen::MatrixXd displacement;
    std::vector<int> permutation;
    // This is potentially a superlattice of the originally passed in parent structure
    xtal::Lattice parent;
    xtal::Lattice child;
};
/// Holds the parameters that are required to conduct a structure map
/// This includes: the parent structure, a point group that determines
/// equivalence,the allowed species on each site of the parent structure,
/// the lattice vs. basis weighting, the maximum allowed volume change from
/// the parent structure, options to the algorithm (sym_basis,sym_strain,robust,strict), tolerance
/// for comparisons, the minimum allowed vacancy concentration, the
/// maximum allowed vacancy concentration, the number of best maps you want,
/// the maximum cost desired, the minimum cost desired, whether or not to keep
/// invalid mappings, whether the structure being mapped is ideal, the
/// potential to impose a lattice to map the test structure onto (must be a superlattice
/// of the parent)
struct MappingInput
{
    typedef CASM::xtal::SymOp SymOp;
    typedef CASM::xtal::SimpleStructure::SpeciesMode SpecMode;
    typedef CASM::xtal::StrucMapping::AllowedSpecies AllowedSpeciesType;

    MappingInput(const casmutils::xtal::Structure& parent)
        : parent(parent),
          mode(SpecMode::ATOM),
          strain_weight(0.5),
          max_volume_change(0.5),
          options(1u << 1),
          tol(CASM::TOL),
          min_va_frac(0.0),
          max_va_frac(0.5),
          num_best_maps(1),
          max_cost(1e20),
          min_cost(-tol),
          keep_invalid_mapping_nodes(false),
          is_ideal(false),
          imposing_lattice(false),
          lattice_to_impose(parent.lattice())
    {
        point_group = CASM::xtal::make_factor_group(parent.__get<CASM::xtal::BasicStructure>());
        for (const auto& site : parent.basis_sites())
        {
            std::vector<std::string> at_site_occs;
            at_site_occs.push_back(site.label());
            allowed_species.push_back(at_site_occs);
        }
    }

    casmutils::xtal::Structure parent;
    std::vector<SymOp> point_group;
    SpecMode mode;
    AllowedSpeciesType allowed_species;
    double strain_weight;
    double max_volume_change;
    int options;
    double tol;
    double min_va_frac;
    double max_va_frac;
    int num_best_maps;
    double max_cost;
    double min_cost;
    bool keep_invalid_mapping_nodes;
    bool is_ideal;
    bool imposing_lattice;
    casmutils::xtal::Lattice lattice_to_impose;
};
/// Can map a structure to its internal reference can be used for mapping many
/// different test structures to the same reference.
class StructureMapper
{
public:
    StructureMapper(const MappingInput& input);

    mapping::MappingNode map(const xtal::Structure& mappable_struc) const;
    mapping::MappingNode ideal_map(const xtal::Structure& mappable_struc) const;

private:
    CASM::xtal::StrucMapper mapper;
};
} // namespace mapping
} // namespace casmutils
#endif
