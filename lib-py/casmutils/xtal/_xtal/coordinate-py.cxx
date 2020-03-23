#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/misc.hpp>
#include <fstream>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace wrappy
{
using namespace casmutils;

namespace Coordinate
{
std::string __str__(const xtal::Coordinate& printable)
{
    std::ostringstream sstream;
    sstream << printable.cart().transpose();
    return sstream.str();
}

/* bool is_equal(const xtal::Coordinate& lhs, const xtal::Coordinate& rhs, double tol) */
/* { */
/*     return casmutils::is_equal<xtal::CoordinateEquals_f>(lhs,rhs,tol); */
/* } */

//TODO: Could have different comparison functions, e.g. periodic compare, and that takes a Lattice
//and a tolerance

} // namespace Coordinate
} // namespace wrappy
