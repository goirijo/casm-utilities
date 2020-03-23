#ifndef LATTICE_PY_HH
#define LATTICE_PY_HH

#include <string>

namespace casmutils::xtal
{
class Lattice;
}

namespace wrappy
{
using namespace casmutils;
namespace Lattice
{
std::string __str__(const xtal::Lattice& printable);
} // namespace Lattice
} // namespace wrappy

#endif
