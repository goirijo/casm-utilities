#include <casmutils/xtal/lattice.hpp>
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
namespace Lattice
{
std::string __str__(const xtal::Lattice& printable)
{
    std::ostringstream sstream;
    sstream << printable.column_vector_matrix();
    return sstream.str();
}
} // namespace Lattice
} // namespace wrappy
