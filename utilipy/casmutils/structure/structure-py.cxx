#include <boost/filesystem.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"
#include "casm/crystallography/Structure.hh"
#include <string>

//******************************************************************************************************//
//******************************************************************************************************//

namespace WrapPy
{
Rewrap::Structure from_poscar(const std::string& filename);
void to_poscar(const Rewrap::Structure& writeable, const std::string& filename);

//******************************************************************************************************//

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

PYBIND11_MODULE(_structure, m)
{
    using namespace pybind11;
    using namespace WrapPy;

    m.doc() = "Raw python bindings for a re-wrapped CASM::Structure class.";

    class_<Rewrap::Structure>(m, "Structure")
        .def("__str__", &__str__)
        .def("is_primitive", &Rewrap::Structure::is_primitive)
        .def("primitive", &Rewrap::Structure::primitive)
        .def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&))Simplicity::make_niggli)
        .def_static("from_poscar", from_poscar)
        .def("to_poscar", to_poscar);

    m.def("make_super_structure", Simplicity::make_super_structure);
    m.def("apply_strain", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::VectorXd&, const std::string&))Simplicity::apply_strain);
    m.def("apply_deformation", (Rewrap::Structure(*)(const Rewrap::Structure&, const Eigen::Matrix3d&))Simplicity::apply_deformation);
    m.def("structure_score", (std::vector<std::pair<double,double>>(*)(const Rewrap::Structure&, const std::vector<Rewrap::Structure>&))Simplicity::structure_score);
}
} // namespace WrapPy
