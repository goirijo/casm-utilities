#include "casmutils/misc.hpp"
#include "casmutils/exceptions.hpp"

namespace extend
{
} // namespace extend

namespace io
{
Eigen::IOFormat coord_format() { return Eigen::IOFormat(4, 0, 0, ", ", "\n", "[", "]"); }
} // namespace io
