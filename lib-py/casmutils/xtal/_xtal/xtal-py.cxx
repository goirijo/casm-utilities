#include "./coordinate-py.hpp"
#include "./lattice-py.hpp"
#include "./rocksalttoggler-py.hpp"
#include "./site-py.hpp"
#include "./structure-py.hpp"
#include "casmutils/xtal/lattice.hpp"

#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/rocksalttoggler.hpp>
#include <casmutils/xtal/structure.hpp>
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
using namespace casmutils;
PYBIND11_MODULE(_xtal, m)
{
    using namespace pybind11;

    m.doc() = "Python bindings for classes and functions related to crystal structures.";

    {
        using namespace wrappy::Structure;
        class_<xtal::Structure>(m, "Structure")
            .def("__str__", __str__)
            //.def("is_primitive", &xtal::Structure::is_primitive)
            .def("make_niggli", (xtal::Structure(*)(const xtal::Structure&))casmutils::xtal::make_niggli)
            .def("set_lattice", set_lattice)
            .def_static("from_poscar", from_poscar)
            .def("to_poscar", to_poscar);
    }

    {
        using namespace wrappy::Lattice;
        class_<xtal::Lattice>(m, "Lattice")
            .def(init<const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&>())
            .def("__str__", __str__)
            .def("a",&xtal::Lattice::a)
            .def("b",&xtal::Lattice::b)
            .def("c",&xtal::Lattice::c);
    }

    {
        typedef xtal::Coordinate xCoord;
        class_<xtal::Coordinate>(m, "Coordinate")
            .def(init<const Eigen::Vector3d&&>())
            .def_static("from_fractional",
                        (xCoord(*)(const Eigen::Vector3d&, const xtal::Lattice&)) & xCoord::from_fractional)
            .def("__str__", wrappy::Coordinate::__str__)
            .def("__add__", &xCoord::operator+, pybind11::is_operator())
            .def("__iadd__", &xCoord::operator+=)
            .def("_bring_within_const", (xCoord(xCoord::*)(const xtal::Lattice&) const) & xCoord::bring_within)
            .def("_bring_within", (void (xCoord::*)(const xtal::Lattice&)) & xCoord::bring_within)
            .def("_cart_const", &xCoord::cart)
            .def("_frac_const", &xCoord::frac);
    }

    {
        class_<xtal::CoordinateEquals_f>(m, "CoordinateEquals_f")
            .def(init<xtal::Coordinate, double>())
            .def("__call__", &xtal::CoordinateEquals_f::operator());
    }

    {
        using namespace wrappy::Site;
        class_<xtal::Site>(m, "Site")
            .def(init<const Eigen::Vector3d&, const std::string&>())
            .def(init<const xtal::Coordinate&, const std::string&>())
            .def("__str__", __str__)
            .def("cart", &xtal::Site::cart)
            .def("frac", &xtal::Site::frac);
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
            .def("activate", (void (RSOT::*)(const xtal::Coordinate&)) & RSOT::activate)
            .def("activate", (void (RSOT::*)(RSOT::index)) & RSOT::activate)
            .def("activate_all", &RSOT::activate_all)
            .def("deactivate", (void (RSOT::*)(const xtal::Coordinate&)) & RSOT::deactivate)
            .def("deactivate", (void (RSOT::*)(RSOT::index)) & RSOT::deactivate)
            .def("deactivate_all", &RSOT::deactivate_all)
            .def("toggle", (void (RSOT::*)(const xtal::Coordinate&)) & RSOT::toggle)
            .def("toggle", (void (RSOT::*)(RSOT::index)) & RSOT::toggle)
            .def("toggle_all", &RSOT::toggle_all)
            .def("nearest_neighbor_distance", &RSOT::nearest_neighbor_distance)
            .def_static("primitive_structure", &RSOT::primitive_structure);
    }

    // clang-format off
    m.def("make_superstructure", casmutils::xtal::make_superstructure);
    m.def("make_primitive", casmutils::xtal::make_primitive);
    m.def("make_niggli", (xtal::Structure(*)(const xtal::Structure&))casmutils::xtal::make_niggli);
    m.def("apply_strain", (xtal::Structure(*)(const xtal::Structure&, const Eigen::VectorXd&, const std::string&))casmutils::xtal::apply_strain);
    m.def("apply_deformation", (xtal::Structure(*)(const xtal::Structure&, const Eigen::Matrix3d&))casmutils::xtal::apply_deformation);
    //m.def("structure_score", (std::vector<std::pair<double, double>>(*)(const rewrap::Structure&, const std::vector<rewrap::Structure>&))casmutils::xtal::structure_score);
    m.def("make_superstructures_of_volume", (std::vector<xtal::Structure>(*)(const xtal::Structure&, const int))casmutils::xtal::make_superstructures_of_volume);
    // clang-format on
}
} // namespace wrappy
