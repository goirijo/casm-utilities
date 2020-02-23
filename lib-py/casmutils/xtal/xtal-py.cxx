#include "casmutils/exceptions.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/xtal/rocksalttoggler.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <string>
#include <fstream>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace wrappy
{
namespace Structure
{
rewrap::Structure from_poscar(const std::string& filename);
void to_poscar(const rewrap::Structure& writeable, const std::string& filename);

std::string __str__(const rewrap::Structure& printable)
{
    std::ostringstream sstream;
    simplicity::print_poscar(printable, sstream);
    return sstream.str();
}

rewrap::Structure from_poscar(const std::string& filename) { return rewrap::Structure::from_poscar(filename); }

void to_poscar(const rewrap::Structure& writeable, const std::string& filename)
{
    simplicity::write_poscar(writeable, filename);
    return;
}

void set_lattice(rewrap::Structure* self, const rewrap::Lattice& new_lattice, std::string mode)
{
    if (mode.size() == 0)
    {
        throw except::BadCoordMode();
    }

    char m = std::tolower(mode[0]);
    switch (m)
    {
    case 'f':
        self->set_lattice(new_lattice, rewrap::COORD_TYPE::FRAC);
        break;

    case 'c':
        self->set_lattice(new_lattice, rewrap::COORD_TYPE::CART);
        break;

    default:
        throw except::BadCoordMode();
    }

    return;
}

} // namespace Structure

namespace Lattice
{
std::string __str__(const rewrap::Lattice& printable)
{
    std::ostringstream sstream;
    sstream << printable.column_vector_matrix();
    return sstream.str();
}
} // namespace Lattice

namespace Site
{
std::string __str__(const ::rewrap::Site& printable)
{
    std::ostringstream sstream;
    throw except::NotImplemented();
    return sstream.str();
}
} // namespace Site

namespace Coordinate
{
std::string __str__(const rewrap::Coordinate& printable)
{
    std::ostringstream sstream;
    sstream << printable.cart().transpose();
    return sstream.str();
}
} // namespace Coordinate

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

PYBIND11_MODULE(_xtal, m)
{
    using namespace pybind11;

    m.doc() = "Python bindings for classes and functions related to crystal structures.";

    {
        using namespace wrappy::Structure;
        class_<rewrap::Structure>(m, "Structure")
            .def("__str__", __str__)
            //.def("is_primitive", &rewrap::Structure::is_primitive)
            .def("make_niggli", (rewrap::Structure(*)(const rewrap::Structure&))simplicity::make_niggli)
            .def("set_lattice", set_lattice)
            .def_static("from_poscar", from_poscar)
            .def("to_poscar", to_poscar);
    }

    {
        using namespace wrappy::Lattice;
        class_<rewrap::Lattice>(m, "Lattice").def(init<const Eigen::Matrix3d&>()).def("__str__", __str__);
    }

    {
        using namespace wrappy::Site;
        class_<rewrap::Site>(m, "Site")
            .def(init<const Eigen::Vector3d&, const std::string&>())
            .def(init<const rewrap::Coordinate&, const std::string&>())
            .def("__str__", __str__)
            .def("cart", &rewrap::Site::cart)
            .def("frac", &rewrap::Site::frac);
    }

    {
        using namespace wrappy::Coordinate;
        class_<rewrap::Coordinate>(m, "Coordinate")
            .def(init<const Eigen::Vector3d&&>())
            .def("__str__", __str__)
            .def("cart", &rewrap::Coordinate::cart)
            .def("frac", &rewrap::Coordinate::frac);
    }

    {
        using namespace wrappy::RockSaltToggler;
        typedef enumeration::RockSaltOctahedraToggler RSOT;
        class_<RSOT>(m, "RockSaltToggler")
            .def("__str__", __str__)
            .def_static("relative_to_primitive", &RSOT::relative_to_primitive)
            .def("all_octahedron_center_coordinates", &RSOT::all_octahedron_center_coordinates)
            .def("to_poscar", to_poscar)
            .def("structure", &RSOT::structure)
            .def("activate", (void (RSOT::*)(const rewrap::Coordinate&)) & RSOT::activate)
            .def("activate", (void (RSOT::*)(RSOT::index)) & RSOT::activate)
            .def("activate_all", &RSOT::activate_all)
            .def("deactivate", (void (RSOT::*)(const rewrap::Coordinate&)) & RSOT::deactivate)
            .def("deactivate", (void (RSOT::*)(RSOT::index)) & RSOT::deactivate)
            .def("deactivate_all", &RSOT::deactivate_all)
            .def("toggle", (void (RSOT::*)(const rewrap::Coordinate&)) & RSOT::toggle)
            .def("toggle", (void (RSOT::*)(RSOT::index)) & RSOT::toggle)
            .def("toggle_all", &RSOT::toggle_all)
            .def("nearest_neighbor_distance", &RSOT::nearest_neighbor_distance)
            .def_static("primitive_structure", &RSOT::primitive_structure);
    }

    // clang-format off
    m.def("make_super_structure", simplicity::make_super_structure);
    m.def("make_primitive", simplicity::make_primitive);
    m.def("make_niggli", (rewrap::Structure(*)(const rewrap::Structure&))simplicity::make_niggli);
    m.def("apply_strain", (rewrap::Structure(*)(const rewrap::Structure&, const Eigen::VectorXd&, const std::string&))simplicity::apply_strain);
    m.def("apply_deformation", (rewrap::Structure(*)(const rewrap::Structure&, const Eigen::Matrix3d&))simplicity::apply_deformation);
    m.def("structure_score", (std::vector<std::pair<double, double>>(*)(const rewrap::Structure&, const std::vector<rewrap::Structure>&))simplicity::structure_score);
    m.def("make_superstructures_of_volume", (std::vector<rewrap::Structure>(*)(const rewrap::Structure&, const int))simplicity::make_superstructures_of_volume);
    m.def("make_boxiest_superstructure_of_volume", (rewrap::Structure(*)(const rewrap::Structure&, const int))simplicity::make_boxiest_superstructure_of_volume);
    // clang-format om
}
} // namespace wrappy
