#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/SimpleStrucMapCalculator.hh>
#include <casm/crystallography/StrucMapping.hh>
#include <casmutils/exceptions.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/sym/cartesian.hpp>

namespace casmutils
{
namespace mapping
{
/// Holds the results of a structure map
/// Fundamentally this includes some strain representation,
/// displacement representation, site assignment matrix/permutation vector,
/// reference superlattice, rigid rotation, and rigid translation.
/// These quantities can be transformed and then converted to a mapping score
/// using strain_cost and basis_cost
struct MappingReport
{
    MappingReport(const CASM::xtal::MappingNode& casm_mapping_node)
        : isometry(casm_mapping_node.isometry()),
          stretch(casm_mapping_node.stretch()),
          translation(casm_mapping_node.translation()),
          displacement(casm_mapping_node.atom_displacement),
          reference_lattice(casm_mapping_node.lat_node.parent.superlattice()),
          mapped_lattice(casm_mapping_node.lat_node.child.superlattice()),
          lattice_cost(casm_mapping_node.lat_node.cost),
          basis_cost(casm_mapping_node.basis_node.cost)
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
    // TODO: !! Explain how the permutation vector works
    std::vector<int> permutation;

    double lattice_cost;
    double basis_cost;

    // This is potentially a superlattice of the originally passed in reference structure
    xtal::Lattice reference_lattice;
    xtal::Lattice mapped_lattice;
};

/// Holds the parameters that are required to conduct a structure map
/// This includes: the reference structure, a point group that determines
/// equivalence,the allowed species on each site of the reference structure,
/// the lattice vs. basis weighting, the maximum allowed volume change from
/// the reference structure, options to the algorithm (sym_basis,sym_strain,robust,strict), tolerance
/// for comparisons, the minimum allowed vacancy concentration, the
/// maximum allowed vacancy concentration, the number of best maps you want,
/// the maximum cost desired, the minimum cost desired, whether or not to keep
/// invalid mappings, whether the structure being mapped is ideal, the
/// potential to impose a lattice to map the test structure onto (must be a superlattice
/// of the reference)
struct MappingInput
{
private:
    typedef CASM::xtal::StrucMapping::AllowedSpecies AllowedSpeciesType;

public:
    // TODO: Rename members with better variables?
    MappingInput(const casmutils::xtal::Structure& reference)
        : reference_structure(reference),
          /* mode(SpecMode::ATOM), */
          strain_weight(0.5),
          max_volume_change(0.5),
          options(1u << 1), //?????
          tol(CASM::TOL),
          min_va_frac(0.0),
          max_va_frac(0.5),
          k_best_maps(1),
          max_cost(1e20),
          min_cost(-tol),
          keep_invalid_mapping_nodes(false),
          impose_reference_lattice(false),
          assume_ideal_lattice(false),
          assume_ideal_structure(false),
          lattice_to_impose(reference.lattice()),
          point_group(CASM::xtal::make_factor_group(reference.__get<CASM::xtal::BasicStructure>()))
    {
        for (const auto& site : reference.basis_sites())
        {
            std::vector<std::string> at_site_occs;
            at_site_occs.push_back(site.label());
            allowed_species.push_back(at_site_occs);
        }
    }

    casmutils::xtal::Structure reference_structure;

    double tol;

    int k_best_maps;
    bool keep_invalid_mapping_nodes;

    double strain_weight;

    double max_volume_change;
    double min_va_frac;
    double max_va_frac;

    double max_cost;
    double min_cost;

    // TODO: Structure point group or lattice point group?
    std::vector<sym::CartOp> point_group;
    bool impose_reference_lattice;
    /// Set to true if the mapped structure has a lattice that is a direct integer tranformation
    /// of the reference, but the basis is relaxed.
    bool assume_ideal_lattice;
    /// Set to true if you know that the mapped structure is a direct integer transformation of the reference.
    /// Implies ideal lattice.
    bool assume_ideal_structure;
    AllowedSpeciesType allowed_species;

    // TODO: Unclear what this could be
    int options;

private:
    // TODO: This might eventually collapse into ATOM mode only, so it's disabled for now
    /* SpecMode mode; */
    casmutils::xtal::Lattice lattice_to_impose;
};

/// Can map a structure to its internal reference can be used for mapping many
/// different test structures to the same reference.
class StructureMapper_f
{
public:
    StructureMapper_f(const MappingInput& input);
    std::vector<MappingReport> operator()(const xtal::Structure& mappable_struc) const;

private:
    CASM::xtal::StrucMapper mapper;
    MappingInput settings;

    std::vector<mapping::MappingReport> map(const xtal::Structure& mappable_struc) const;
    std::vector<mapping::MappingReport> ideal_map(const xtal::Structure& mappable_struc) const;
};

/// Calculates lattice and basis score from ideal lattice, stretch tensor and displacement matrix
/// Returns scores for lattice (first) and basis (second) as a pair.
std::pair<double, double> structure_score(const mapping::MappingReport& mapping_data);

/// Map a single structure onto a reference structure with default settings
std::vector<mapping::MappingReport> map_structure(const xtal::Structure& map_reference_struc, const xtal::Structure& mappable_struc);

} // namespace mapping
} // namespace casmutils
#endif
