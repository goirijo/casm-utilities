#ifndef SIMPLICITY_HH
#define SIMPLICITY_HH

#include <casm/CASM_global_definitions.hh>
#include <iosfwd>

namespace CASM
{
class Structure;
class PrimClex;
class Site;

template <class T>
class BasicStructure;
}

// Functions to make dealing with the casm interface less of a PITA
namespace simple
{
/// Go through the motions of creating the VaspIO object and print a normal POSCAR format
void structure_print(std::ostream &stream, const CASM::Structure &struc);

/// Given a path to prim.json like file, create an empty PrimClex
CASM::PrimClex primclex_from_path(CASM::fs::path prim_path);

/// Given a path to a POS.vasp like file, create a BasicStructure
CASM::BasicStructure<CASM::Site> basic_structure_from_path(CASM::fs::path struc_path);

/// Use zero intelligence defaults and discard extra information to import BasicStructure into PrimClex
bool primclex_import_basic_structure(CASM::PrimClex &pclex, CASM::BasicStructure<CASM::Site> new_struc);

/// Round a MatrixXd to a MatrixXi without worrying about casting
Eigen::MatrixXi easy_matrix_round(const Eigen::MatrixXd &mat);
}

#endif
