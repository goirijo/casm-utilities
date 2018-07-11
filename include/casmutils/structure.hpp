#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include <iostream>
#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>

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

    Structure(CASM::Structure init_struc);

    bool is_primitive() const;

    Structure primitive() const;

private:
};
}

namespace Simplicity
{
void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename);
void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream);
Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli);
Rewrap::Structure make_niggli(Rewrap::Structure non_niggli);
void make_niggli(Rewrap::Structure* non_niggli);
Rewrap::Structure make_primitive(const Rewrap::Structure& input);
}


#endif
