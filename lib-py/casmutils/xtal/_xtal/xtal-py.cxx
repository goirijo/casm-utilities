#include "./coordinate-py.hpp"
#include "./lattice-py.hpp"
#include "./rocksalttoggler-py.hpp"
#include "./site-py.hpp"
#include "./structure-py.hpp"

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
