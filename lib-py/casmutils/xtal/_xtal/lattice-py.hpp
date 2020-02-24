#ifndef LATTICE_PY_HH
#define LATTICE_PY_HH

#include <string>

namespace rewrap
{
    class Lattice;
}

namespace wrappy
{
namespace Lattice
{
std::string __str__(const rewrap::Lattice& printable);
} // namespace Lattice
} // namespace wrappy


#endif
