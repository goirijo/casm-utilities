#include "./coordinate-py.hpp"
#include "./lattice-py.hpp"
#include "./rocksalttoggler-py.hpp"
#include "./site-py.hpp"
#include "./structure-py.hpp"
#include "casmutils/xtal/lattice.hpp"

#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/rocksalttoggler.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/symmetry.hpp>
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
            .def(init<const xtal::Lattice&, const std::vector<xtal::Site>&>())
            .def("__str__", __str__)
            .def_static("_from_poscar", from_poscar)
            .def("_to_poscar", to_poscar)
            .def("_lattice_const", &xtal::Structure::lattice)
            .def("_set_lattice_const", set_lattice_const)
            .def("_set_lattice", set_lattice)
            .def("_within", &xtal::Structure::within)
            .def("_basis_sites_const", &xtal::Structure::basis_sites)
            .def("make_niggli", pybind11::overload_cast<const xtal::Structure&>(casmutils::xtal::make_niggli));
    }

    {
        using namespace wrappy::Lattice;
        class_<xtal::Lattice>(m, "Lattice")
            .def(init<const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&>())
            .def(init<const Eigen::Matrix3d&>())
            .def("__str__", __str__)
            .def("a", &xtal::Lattice::a)
            .def("b", &xtal::Lattice::b)
            .def("c", &xtal::Lattice::c)
            .def("__getitem__", &xtal::Lattice::operator[])
            .def("volume", &xtal::Lattice::volume)
            .def("column_vector_matrix", &xtal::Lattice::column_vector_matrix)
            .def("row_vector_matrix", &xtal::Lattice::row_vector_matrix)
            .def("make_niggli", pybind11::overload_cast<const xtal::Lattice&>(casmutils::xtal::make_niggli));
    }

    {
        class_<xtal::LatticeEquals_f>(m, "LatticeEquals_f")
            .def(init<xtal::Lattice, double>())
            .def("__call__", &xtal::LatticeEquals_f::operator());
    }

    {
        using namespace wrappy::Coordinate;
        class_<xtal::Coordinate>(m, "Coordinate")
            .def(init<const Eigen::Vector3d&&>())
            .def_static("from_fractional",
                        pybind11::overload_cast<const Eigen::Vector3d&, const xtal::Lattice&>(
                            &xtal::Coordinate::from_fractional))
            .def("__str__", __str__)
            .def("__add__", &xtal::Coordinate::operator+, pybind11::is_operator())
            .def("__iadd__", &xtal::Coordinate::operator+=)
            .def("_bring_within_const",
                 pybind11::overload_cast<const xtal::Lattice&>(&xtal::Coordinate::bring_within, pybind11::const_))
            .def("_bring_within", pybind11::overload_cast<const xtal::Lattice&>(&xtal::Coordinate::bring_within))
            .def("_cart_const", &xtal::Coordinate::cart)
            .def("_frac_const", &xtal::Coordinate::frac)
            .def("__rmul__", [](const xtal::Coordinate& coord, const sym::CartOp& symop) {
                return symop * coord;
            });
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
            .def("_cart_const", &xtal::Site::cart)
            .def("_frac_const", &xtal::Site::frac)
            .def("_label_const", &xtal::Site::label)
            .def("__rmul__", [](const xtal::Site& site, const sym::CartOp& symop) {
                return symop * site;
            });
    }

    {
        class_<xtal::SiteEquals_f>(m, "SiteEquals_f")
            .def(init<xtal::Site, double>())
            .def("__call__", &xtal::SiteEquals_f::operator());
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
    m.def("make_niggli", (xtal::Lattice(*)(const xtal::Lattice&))casmutils::xtal::make_niggli);
    m.def("apply_strain", (xtal::Structure(*)(const xtal::Structure&, const Eigen::VectorXd&, const std::string&))casmutils::xtal::apply_strain);
    m.def("apply_deformation", (xtal::Structure(*)(const xtal::Structure&, const Eigen::Matrix3d&))casmutils::xtal::apply_deformation);
    m.def("make_superstructures_of_volume", (std::vector<xtal::Structure>(*)(const xtal::Structure&, const int))casmutils::xtal::make_superstructures_of_volume);

    m.def("make_point_group", casmutils::xtal::make_point_group);
    m.def("make_factor_group", casmutils::xtal::make_factor_group);
	m.def("_symmetrize_lattice",(xtal::Lattice(*)(const xtal::Lattice&, const std::vector<sym::CartOp>&))casmutils::xtal::symmetrize);
	m.def("_symmetrize_structure",(xtal::Structure(*)(const xtal::Structure&, const std::vector<sym::CartOp>&))casmutils::xtal::symmetrize);
    // clang-format on
}
} // namespace wrappy
