#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>

namespace CASM
{
class Structure;
}

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
void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename);
void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream);
Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli);
Rewrap::Structure make_niggli(Rewrap::Structure non_niggli);
void make_niggli(Rewrap::Structure* non_niggli);
Rewrap::Structure make_primitive(const Rewrap::Structure& input);
} // namespace Simplicity

#endif
