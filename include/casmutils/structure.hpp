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
} // namespace Simplicity

#endif
