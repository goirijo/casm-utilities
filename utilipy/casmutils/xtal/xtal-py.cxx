#include <boost/filesystem.hpp>
#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/lattice.hpp"
#include "casmutils/structure_tools.hpp"
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
std::string __str__(const Rewrap::Site& printable)
{
    std::ostringstream sstream;
    sstream << printable.cart().transpose()<<"    "<<printable.current_occupant_name();
    return sstream.str();
}
} // namespace Site

PYBIND11_MODULE(_xtal, m)
{
    using namespace pybind11;

    m.doc() = "Raw python bindings for a re-wrapped CASM::Structure class.";

    {
        using namespace WrapPy::Structure;
        class_<Rewrap::Structure>(m, "Structure")
            .def("__str__", __str__)
            .def("is_primitive", &Rewrap::Structure::is_primitive)
            .def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&))Simplicity::make_niggli)
            .def_static("from_poscar", from_poscar)
            .def("to_poscar", to_poscar);
    }

    {
        using namespace WrapPy::Lattice;
        class_<Rewrap::Lattice>(m, "Lattice")
            .def(init<const Eigen::Matrix3d&>())
            .def("__str__", __str__);
    }

    {
        using namespace WrapPy::Site;
        class_<Rewrap::Site>(m, "Site")
            .def(init<const Eigen::Vector3d&, const std::vector<std::string>&>())
            .def("__str__", __str__)
            .def("cart", &Rewrap::Site::cart)
            .def("frac", &Rewrap::Site::frac)
            .def("current_occupant_name", &Rewrap::Site::current_occupant_name);
    }

    m.def("make_super_structure", Simplicity::make_super_structure);
    m.def("make_primitive", Simplicity::make_primitive);
    m.def("apply_strain", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::VectorXd&, const std::string&))Simplicity::apply_strain);
    m.def("apply_deformation", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::Matrix3d&))Simplicity::apply_deformation);
    m.def("structure_score", (std::vector<std::pair<double, double>>(*)(const Rewrap::Structure&, const std::vector<Rewrap::Structure>&))Simplicity::structure_score);
    m.def("make_superstructures_of_volume", (std::vector<Rewrap::Structure>(*)(const Rewrap::Structure&, const int))Simplicity::make_superstructures_of_volume);
    m.def("make_boxiest_superstructure_of_volume", (Rewrap::Structure(*)(const Rewrap::Structure&, const int))Simplicity::make_boxiest_superstructure_of_volume);
}
} // namespace WrapPy
