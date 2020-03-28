#ifndef STRUCTURE_TOOLS_HH
#define STRUCTURE_TOOLS_HH

#include <casmutils/xtal/structure.hpp>

namespace casmutils
{

namespace mapping
{
struct MappingReport;
struct MappingInput;
} // namespace mapping

namespace xtal
{ /// Given a Structure, write out its information into a file in a vasp compatible format
void write_poscar(const Structure& printable, const fs::path& filename);

/// Given a Structure, print out its information to the given stream in a vasp compatible format
void print_poscar(const Structure& printable, std::ostream& outstream);

/// Return a copy of the given Structure that has been converted to its standard niggli form
Structure make_niggli(const Structure& non_niggli);

/// Modify the given Structure to standard niggli form
void make_niggli(Structure* non_niggli);

/// Return a Structure that is the primitive of the provided one
Structure make_primitive(const Structure& input);

/// Returns a super structure after applying a transformation matrix to the structure.
/// Applies transformation to the lattice and uses CASM::Structure::create_superstruc to fill the basis.
/// transformed_lattice =  original_lattice * transformation_matrix
Structure make_superstructure(const Structure& struc, const Eigen::Matrix3i& col_transf_mat);

/// Takes a pointer to a structure and applies deformation to that structure.
/// Deforms the lattice and keeps the basis constant in fractional space.
/// deformed_lattice = deformation_tensor * undeformed_lattice.
/// Does not remove rotation even if it exists!!!
void apply_deformation(Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor);
Structure apply_deformation(const Structure& struc_ptr, const Eigen::Matrix3d& deformation_tensor);

/// Takes the unique elements of a strain metric and arranges them correctly
/// in a matrix with the proper normalization
Eigen::Matrix3d rollup_strain_metric(const Eigen::Ref<const Eigen::VectorXd>& unrolled_strain);

/// Takes a pointer to a structure and applies strain to that structure.
/// Input is unrolled strain in conventional metrics as defined in the mode.
/// Allowed modes are 'GL' [Green-Lagrange], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky] and throws an error if
/// the mode is not in this list. Uses functions from CASM::StrainConverter class to roll up the strain and obtain a
/// deformation tensor. Applies deformation using apply_deformation function.
void apply_strain(Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode);
Structure apply_strain(const Structure& struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode);

/// Given a structure, find all the superstructures between volumes min_vol and max_vol
std::vector<Structure> make_superstructures_of_volume(const Structure& structure, const int volume);

///// Map a vector of structures onto a single reference structures, return a vector of score pairs
///// for the lattice (first) and basis (second).
// std::vector<std::pair<double, double>> structure_score(const Structure& map_reference_struc,
//                                                       const std::vector<Structure>& mappable_struc_vec);

/// Calculates lattice and basis score from ideal lattice, stretch tensor and displacement matrix
/// Returns scores for lattice (first) and basis (second) as a pair.
std::pair<double, double> structure_score(const mapping::MappingReport& mapping_data);

/// Map a single structure onto a reference structure with default settings
std::vector<mapping::MappingReport> map_structure(const Structure& map_reference_struc, const Structure& mappable_struc);

} // namespace xtal
} // namespace casmutils

#endif
