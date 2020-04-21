#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
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
namespace Structure
{
std::string __str__(const xtal::Structure& printable)
{
    std::ostringstream sstream;
    casmutils::xtal::print_poscar(printable, sstream);
    return sstream.str();
}

xtal::Structure from_poscar(const std::string& filename) { return xtal::Structure::from_poscar(filename); }

void to_poscar(const xtal::Structure& writeable, const std::string& filename)
{
    casmutils::xtal::write_poscar(writeable, filename);
    return;
}

xtal::COORD_TYPE string_to_coord_mode(std::string mode)
{
    if (mode.size() == 0)
    {
        throw except::BadCoordMode();
    }

    char m = std::tolower(mode[0]);

    switch (m)
    {

    case 'f':
        return xtal::COORD_TYPE::FRAC;

    case 'c':
        return xtal::COORD_TYPE::CART;

    default:
        throw except::BadCoordMode();
    }
}

void set_lattice(xtal::Structure* self, const xtal::Lattice& new_lattice, std::string mode)
{
    self->set_lattice(new_lattice, string_to_coord_mode(mode));
    return;
}

xtal::Structure set_lattice_const(xtal::Structure* self, const xtal::Lattice& new_lattice, std::string mode)
{
    const xtal::Structure self_structure = *self;
    return self_structure.set_lattice(new_lattice, string_to_coord_mode(mode));
}
} // namespace Structure

} // namespace wrappy
