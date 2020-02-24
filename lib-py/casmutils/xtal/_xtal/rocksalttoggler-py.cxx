#include "./structure-py.hpp"
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/rocksalttoggler.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <fstream>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace wrappy
{
namespace RockSaltToggler
{
void to_poscar(const enumeration::RockSaltOctahedraToggler& writeable, const std::string& filename)
{
    std::ofstream rs_outstream;
    rs_outstream.open(filename);
    writeable.print(rs_outstream);
    rs_outstream.close();

    return;
}

std::string __str__(const enumeration::RockSaltOctahedraToggler& printable)
{
    auto structure = printable.structure();
    return Structure::__str__(structure);
}
} // namespace RockSaltToggler

} // namespace wrappy
