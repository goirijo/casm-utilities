#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include <iostream>

namespace CASM
{
    namespace fs
    {
        class path;
    }

    class Structure;
}

namespace Simplicity
{
void print_poscar(const CASM::Structure& printable, const CASM::fs::path& filename);
void print_poscar(const CASM::Structure& printable, std::ostream& outstream);
CASM::Structure make_niggli(const CASM::Structure& non_niggli);
CASM::Structure make_primitive(const CASM::Structure& input);
}


#endif
