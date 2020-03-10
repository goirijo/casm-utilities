#include <casmutils/exceptions.hpp>
#include <casmutils/xtal/site.hpp>
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
namespace Site
{
std::string __str__(const xtal::Site& printable)
{
    std::ostringstream sstream;
    throw except::NotImplemented();
    return sstream.str();
}
} // namespace Site
} // namespace wrappy
