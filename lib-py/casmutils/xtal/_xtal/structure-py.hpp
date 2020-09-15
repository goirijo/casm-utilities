#ifndef STRUCTURE_PY_HH
#define STRUCTURE_PY_HH

#include <string>

namespace casmutils
{
namespace xtal
{
class Structure;
class Lattice;
} // namespace xtal
} // namespace casmutils

namespace wrappy
{
using namespace casmutils;
namespace Structure
{
xtal::Structure from_poscar(const std::string& filename);

void to_poscar(const xtal::Structure& writeable, const std::string& filename);

std::string __str__(const xtal::Structure& printable);

void set_lattice(xtal::Structure* self, const xtal::Lattice& new_lattice, std::string mode);

xtal::Structure set_lattice_const(xtal::Structure* self, const xtal::Lattice& new_lattice, std::string mode);

} // namespace Structure

} // namespace wrappy

#endif
