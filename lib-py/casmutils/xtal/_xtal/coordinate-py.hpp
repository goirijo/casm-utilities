#ifndef COORDINATE_PY_HH
#define COORDINATE_PY_HH

#include <string>

namespace casmutils::xtal
{
class Coordinate;
}

namespace wrappy
{
using namespace casmutils;
namespace Coordinate
{
std::string __str__(const xtal::Coordinate& printable);
}
} // namespace wrappy

#endif
