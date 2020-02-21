#include "casmutils/exceptions.hpp"
#include "casmutils/lattice.hpp"
#include "casmutils/misc.hpp"
#include "casmutils/rocksalttoggler.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/structure_tools.hpp"
#include <boost/filesystem.hpp>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace WrapPy
{
namespace Structure
{
Rewrap::Structure from_poscar(const std::string& filename);
void to_poscar(const Rewrap::Structure& writeable, const std::string& filename);

std::string __str__(const Rewrap::Structure& printable)
{
    std::ostringstream sstream;
    Simplicity::print_poscar(printable, sstream);
    return sstream.str();
}

Rewrap::Structure from_poscar(const std::string& filename) { return Rewrap::Structure::from_poscar(filename); }

void to_poscar(const Rewrap::Structure& writeable, const std::string& filename)
{
    Simplicity::write_poscar(writeable, filename);
    return;
}

void set_lattice(Rewrap::Structure* self, const Rewrap::Lattice& new_lattice, std::string mode)
{
    if (mode.size() == 0)
    {
        throw UtilExcept::BadCoordMode();
    }

    char m = std::tolower(mode[0]);
    switch (m)
    {
    case 'f':
        self->set_lattice(new_lattice, Rewrap::COORD_TYPE::FRAC);
        break;

    case 'c':
        self->set_lattice(new_lattice, Rewrap::COORD_TYPE::CART);
        break;

    default:
        throw UtilExcept::BadCoordMode();
    }

    return;
}

} // namespace Structure

namespace Lattice
{
std::string __str__(const Rewrap::Lattice& printable)
{
    std::ostringstream sstream;
    sstream << printable.lat_column_mat();
    return sstream.str();
}
} // namespace Lattice

namespace Site
{
std::string __str__(const ::Rewrap::Site& printable)
{
    std::ostringstream sstream;
    throw UtilExcept::NotImplemented();
    return sstream.str();
}
} // namespace Site

namespace Coordinate
{
std::string __str__(const Rewrap::Coordinate& printable)
{
    std::ostringstream sstream;
    sstream << printable.cart().transpose();
    return sstream.str();
}
} // namespace Coordinate

namespace RockSaltToggler
{
void to_poscar(const SpecializedEnumeration::RockSaltOctahedraToggler& writeable, const std::string& filename)
{
    std::ofstream rs_outstream;
    rs_outstream.open(filename);
    writeable.print(rs_outstream);
    rs_outstream.close();

    return;
}

std::string __str__(const SpecializedEnumeration::RockSaltOctahedraToggler& printable)
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
        using namespace WrapPy::Structure;
        class_<Rewrap::Structure>(m, "Structure")
            .def("__str__", __str__)
            .def("is_primitive", &Rewrap::Structure::is_primitive)
            .def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&))Simplicity::make_niggli)
            .def("set_lattice", set_lattice)
            .def_static("from_poscar", from_poscar)
            .def("to_poscar", to_poscar);
    }

    {
        using namespace WrapPy::Lattice;
        class_<Rewrap::Lattice>(m, "Lattice").def(init<const Eigen::Matrix3d&>()).def("__str__", __str__);
    }

    {
        using namespace WrapPy::Site;
        class_<Rewrap::Site>(m, "Site")
            .def(init<const Eigen::Vector3d&, const std::vector<std::string>&>())
            .def(init<const Rewrap::Coordinate&, const std::vector<std::string>&>())
            .def("__str__", __str__)
            .def("cart", &Rewrap::Site::cart)
            .def("frac", &Rewrap::Site::frac);
    }

    {
        using namespace WrapPy::Coordinate;
        class_<Rewrap::Coordinate>(m, "Coordinate")
            .def(init<const Eigen::Vector3d&&>())
            .def("__str__", __str__)
            .def("cart", &Rewrap::Coordinate::cart)
            .def("frac", &Rewrap::Coordinate::frac);
    }

    {
        using namespace WrapPy::RockSaltToggler;
        typedef SpecializedEnumeration::RockSaltOctahedraToggler RSOT;
        class_<RSOT>(m, "RockSaltToggler")
            .def("__str__", __str__)
            .def_static("relative_to_primitive", &RSOT::relative_to_primitive)
            .def("all_octahedron_center_coordinates", &RSOT::all_octahedron_center_coordinates)
            .def("to_poscar", to_poscar)
            .def("structure", &RSOT::structure)
            .def("activate", (void (RSOT::*)(const Rewrap::Coordinate&)) & RSOT::activate)
            .def("activate", (void (RSOT::*)(RSOT::index)) & RSOT::activate)
            .def("activate_all", &RSOT::activate_all)
            .def("deactivate", (void (RSOT::*)(const Rewrap::Coordinate&)) & RSOT::deactivate)
            .def("deactivate", (void (RSOT::*)(RSOT::index)) & RSOT::deactivate)
            .def("deactivate_all", &RSOT::deactivate_all)
            .def("toggle", (void (RSOT::*)(const Rewrap::Coordinate&)) & RSOT::toggle)
            .def("toggle", (void (RSOT::*)(RSOT::index)) & RSOT::toggle)
            .def("toggle_all", &RSOT::toggle_all)
            .def("nearest_neighbor_distance", &RSOT::nearest_neighbor_distance)
            .def_static("primitive_structure", &RSOT::primitive_structure);
    }

    // clang-format off
    m.def("make_super_structure", Simplicity::make_super_structure);
    m.def("make_primitive", Simplicity::make_primitive);
    m.def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&))Simplicity::make_niggli);
    m.def("apply_strain", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::VectorXd&, const std::string&))Simplicity::apply_strain);
    m.def("apply_deformation", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::Matrix3d&))Simplicity::apply_deformation);
    m.def("structure_score", (std::vector<std::pair<double, double>>(*)(const Rewrap::Structure&, const std::vector<Rewrap::Structure>&))Simplicity::structure_score);
    m.def("make_superstructures_of_volume", (std::vector<Rewrap::Structure>(*)(const Rewrap::Structure&, const int))Simplicity::make_superstructures_of_volume);
    m.def("make_boxiest_superstructure_of_volume", (Rewrap::Structure(*)(const Rewrap::Structure&, const int))Simplicity::make_boxiest_superstructure_of_volume);
    // clang-format om
}
} // namespace WrapPy
