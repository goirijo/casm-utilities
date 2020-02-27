#ifndef ROCKSALTTOGGLER_PY_HH
#define ROCKSALTTOGGLER_PY_HH

#include <string>

namespace enumeration
{
class RockSaltOctahedraToggler;
}

namespace wrappy
{
namespace RockSaltToggler
{
void to_poscar(const enumeration::RockSaltOctahedraToggler& writeable, const std::string& filename);

std::string __str__(const enumeration::RockSaltOctahedraToggler& printable);
} // namespace RockSaltToggler
} // namespace wrappy
#endif
