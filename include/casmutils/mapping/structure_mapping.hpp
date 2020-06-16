#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/SimpleStrucMapCalculator.hh>
#include <casm/crystallography/StrucMapping.hh>
#include <casmutils/exceptions.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>

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
          reference_lattice(casm_mapping_node.lattice_node.parent.superlattice()),
          mapped_lattice(casm_mapping_node.lattice_node.child.superlattice()),
          lattice_cost(casm_mapping_node.lattice_node.cost),
          basis_cost(casm_mapping_node.atomic_node.cost),
          cost(casm_mapping_node.cost)
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
    double cost;

    // This is potentially a superlattice of the originally passed in reference structure
    xtal::Lattice reference_lattice;
    xtal::Lattice mapped_lattice;
};

/// Holds the parameters that are required to conduct a structure map, including
/// the lattice vs. basis weighting, the maximum allowed volume change from
/// the reference structure, options to the algorithm (sym_basis,sym_strain,robust,strict), tolerance
/// for comparisons, the minimum allowed vacancy concentration, the
/// maximum allowed vacancy concentration, the number of best maps you want,
/// the maximum cost desired, the minimum cost desired, whether or not to keep
/// invalid mappings, whether the structure being mapped is ideal, the
/// potential to impose a lattice to map the test structure onto (must be a superlattice
/// of the reference)
// TODO: Explain each parameter in detail
struct MappingInput
{
public:
    MappingInput()
        : /* mode(SpecMode::ATOM), */
          strain_weight(0.5),
          max_volume_change(0.5),
          options(CASM::xtal::StrucMapper::Options::robust),
          tol(CASM::TOL),
          min_vacancy_fraction(0.0),
          max_vacancy_fraction(0.5),
          k_best_maps(1),
          max_cost(1e20),
          min_cost(-tol),
          keep_invalid_mapping_nodes(false),
          impose_reference_lattice(false),
          assume_ideal_lattice(false),
          assume_ideal_structure(false),
          /* assume_deformed_structure(false), */
          use_crystal_symmetry(false)
    {
    }

    double tol;

    int k_best_maps;
    bool keep_invalid_mapping_nodes;

    double strain_weight;

    double max_volume_change;
    double min_vacancy_fraction;
    double max_vacancy_fraction;

    double max_cost;
    double min_cost;

    bool impose_reference_lattice;
    /// Set to true if the mapped structure has a lattice that is a direct integer tranformation
    /// of the reference, but the basis is relaxed.
    bool assume_ideal_lattice;
    /// Set to true if you know that the mapped structure is a direct integer transformation of the reference.
    /// Implies ideal lattice.
    bool assume_ideal_structure;
    /// TODO: More docs
    /* bool assume_deformed_structure; */

    // TODO: Unclear what this could be
    int options;

    /// When true, the factor group of the reference structure is applied to the mapped structure
    /// when performing the mapping
    bool use_crystal_symmetry;

private:
    // TODO: This might eventually collapse into ATOM mode only, so it's disabled for now
    /* SpecMode mode; */
};

/// Can map a structure to its internal reference can be used for mapping many
/// different test structures to the same reference.
/// Default values for the point group is the factor group of the reference structure,
/// and the allowed species are whatever is residing at the reference structure.
// TODO: Explain each constructor argument in detail
class StructureMapper_f
{
public:
    typedef CASM::xtal::StrucMapping::AllowedSpecies AllowedSpeciesType;

    StructureMapper_f(const xtal::Structure& reference,
                      const MappingInput& input,
                      const std::vector<sym::CartOp>& factor_group = {},
                      const AllowedSpeciesType& allowed_species = {});

    // Having this allows passing either factor_group OR allowed_species, both, or neither
    StructureMapper_f(const xtal::Structure& reference,
                      const MappingInput& input,
                      const AllowedSpeciesType& allowed_species)
        : StructureMapper_f(reference, input, {}, allowed_species)
    {
    }

    std::vector<MappingReport> operator()(const xtal::Structure& mappable_struc) const;

private:
    xtal::Structure reference_structure;
    xtal::Lattice lattice_to_impose;
    MappingInput settings;

    std::vector<sym::CartOp> factor_group;
    AllowedSpeciesType allowed_species;

    CASM::xtal::StrucMapper mapper;

    std::vector<mapping::MappingReport> map(const xtal::Structure& mappable_struc) const;
    std::vector<mapping::MappingReport> ideal_map(const xtal::Structure& mappable_struc) const;

    /// Returns the factor group of the reference structure
    std::vector<sym::CartOp> make_default_factor_group() const;

    /// Returns the current species of the reference structure
    AllowedSpeciesType make_default_allowed_species() const;
};

/// Calculates lattice and basis score from ideal lattice, stretch tensor and displacement matrix
/// Returns scores for lattice (first) and basis (second) as a pair.
std::pair<double, double> structure_score(const mapping::MappingReport& mapping_data);

/// Map a single structure onto a reference structure with default settings
std::vector<mapping::MappingReport> map_structure(const xtal::Structure& map_reference_struc,
                                                  const xtal::Structure& mappable_struc);

/// Calculates the symmetry preserving part of the MappingReport according to the given factor group
mapping::MappingReport symmetry_preserving_mapping_report(const mapping::MappingReport& mapping_data,
                                                          const std::vector<sym::CartOp>& group_as_operations,
                                                          const std::vector<sym::PermRep>& group_as_permutations);

} // namespace mapping
} // namespace casmutils
#endif
