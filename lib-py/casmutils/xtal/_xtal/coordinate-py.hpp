#ifndef COORDINATE_PY_HH
#define COORDINATE_PY_HH

#include <string>

namespace rewrap
{
class Coordinate;
}

namespace wrappy
{
namespace Coordinate
{
std::string __str__(const rewrap::Coordinate& printable);
}
} // namespace wrappy

#endif
