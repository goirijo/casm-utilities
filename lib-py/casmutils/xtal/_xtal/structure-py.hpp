#ifndef STRUCTURE_PY_HH
#define STRUCTURE_PY_HH

#include <string>

namespace rewrap
{
    class Structure;
    class Lattice;
}

namespace wrappy
{
namespace Structure
{
rewrap::Structure from_poscar(const std::string& filename);

void to_poscar(const rewrap::Structure& writeable, const std::string& filename);

std::string __str__(const rewrap::Structure& printable);

void set_lattice(rewrap::Structure* self, const rewrap::Lattice& new_lattice, std::string mode);

} // namespace Structure

} // namespace wrappy

#endif
