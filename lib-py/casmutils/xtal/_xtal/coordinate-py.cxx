#include <casmutils/xtal/coordinate.hpp>
#include <fstream>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace wrappy
{

namespace Coordinate
{
std::string __str__(const rewrap::Coordinate& printable)
{
    std::ostringstream sstream;
    sstream << printable.cart().transpose();
    return sstream.str();
}
} // namespace Coordinate
} // namespace wrappy
