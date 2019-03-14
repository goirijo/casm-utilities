#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>

namespace Rewrap
{

class Structure : public CASM::Structure
{
public:
    Structure() = delete;

    /// Construct by providing a path to a POSCAR like file
    static Structure from_poscar(const fs::path& poscar_path);

    /// Construct from parent CASM class
    Structure(CASM::Structure init_struc);
    Structure(Rewrap::fs::path& filename);

    /// Returns true if the structure is already primitive
    bool is_primitive() const;

    /// Creates a new structure that is the primitive cell of *this
    Structure primitive() const;

private:
};
} // namespace Rewrap

namespace Simplicity
{
/// Given a Structure, write out its information into a file in a vasp compatible format
void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename);

/// Given a Structure, print out its information to the given stream in a vasp compatible format
void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream);

/// Return a copy of the given Structure that has been converted to its standard niggli form
Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli);

/// Modify the given Structure to standard niggli form
void make_niggli(Rewrap::Structure* non_niggli);

/// Return a Structure that is the primitive of the provided one
Rewrap::Structure make_primitive(const Rewrap::Structure& input);

/// Returns a super structure after applying a transformation matrix to the structure.
/// Applies transformation to the lattice and uses CASM::Structure::create_superstruc to fill the basis.
/// transformed_lattice =  original_lattice * transformation_matrix
Rewrap::Structure make_super_structure(const Rewrap::Structure& struc, const Eigen::Matrix3i& col_transf_mat);

/// Takes a pointer to a structure and applies deformation to that structure.
/// Deforms the lattice and keeps the basis constant in fractional space.
/// deformed_lattice = deformation_tensor * undeformed_lattice.
/// Does not remove rotation even if it exists!!!
void apply_deformation(Rewrap::Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor);
Rewrap::Structure apply_deformation(const Rewrap::Structure& struc_ptr, const Eigen::Matrix3d& deformation_tensor);

/// Takes a pointer to a structure and applies strain to that structure.
/// Input is unrolled strain in conventional metrics as defined in the mode.
/// Allowed modes are 'GL' [Green-Lagrange], 'EA' [Euler-Almansi], 'B' [Biot], or 'H' [Hencky] and throws an error if
/// the mode is not in this list. Uses functions from CASM::StrainConverter class to roll up the strain and obtain a
/// deformation tensor. Applies deformation using apply_deformation function.
void apply_strain(Rewrap::Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode);
Rewrap::Structure apply_strain(const Rewrap::Structure& struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode);

/// Map a vector of structures onto a single reference structures, return a vector of score pairs
/// for the lattice (first) and basis (second).
std::vector<std::pair<double,double>> structure_score(const Rewrap::Structure& map_reference_struc,
                       const std::vector<Rewrap::Structure>& mappable_struc_vec);

/// Map a single structure onto a reference structure.
/// Returns scores for lattice (first) and basis (second) as a pair.
std::pair<double,double> structure_score(const Rewrap::Structure& map_reference_struc,
                       const Rewrap::Structure& mappable_struc);

} // namespace Simplicity

#endif
