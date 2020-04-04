#include <casmutils/mapping/structure_mapping.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <fstream>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//******************************************************************************************************//
//******************************************************************************************************//

namespace casmutils
{
}

namespace wrappy
{

using namespace casmutils;
PYBIND11_MODULE(_mapping, m)
{
    using namespace pybind11;

    m.doc() = "Python bindings for classes and functions related to mapping objects from one to another\
               , such as crystal structures or clusters.";

    {
        class_<mapping::MappingReport>(m, "MappingReport")
            .def_readonly("isometry",&mapping::MappingReport::isometry)
            .def_readonly("stretch",&mapping::MappingReport::stretch)
            .def_readonly("translation",&mapping::MappingReport::translation)
            .def_readonly("displacement",&mapping::MappingReport::displacement)
            .def_readonly("permutation",&mapping::MappingReport::permutation)
            .def_readonly("lattice_cost",&mapping::MappingReport::lattice_cost)
            .def_readonly("basis_cost",&mapping::MappingReport::basis_cost)
            .def_readonly("reference_lattice",&mapping::MappingReport::reference_lattice)
            .def_readonly("mapped_lattice",&mapping::MappingReport::mapped_lattice);
    }

    {
        class_<mapping::MappingInput>(m, "MappingInput")
            .def(init<>())
    .def_readwrite("tol",&mapping::MappingInput::tol)
    .def_readwrite("k_best_maps",&mapping::MappingInput::k_best_maps)
    .def_readwrite("keep_invalid_mapping_nodes",&mapping::MappingInput::keep_invalid_mapping_nodes)
    .def_readwrite("strain_weight",&mapping::MappingInput::strain_weight)
    .def_readwrite("max_volume_change",&mapping::MappingInput::max_volume_change)
    .def_readwrite("min_va_frac",&mapping::MappingInput::min_va_frac)
    .def_readwrite("max_va_frac",&mapping::MappingInput::max_va_frac)
    .def_readwrite("max_cost",&mapping::MappingInput::max_cost)
    .def_readwrite("min_cost",&mapping::MappingInput::min_cost)
    .def_readwrite("impose_reference_lattice",&mapping::MappingInput::impose_reference_lattice)
    .def_readwrite("assume_ideal_lattice",&mapping::MappingInput::assume_ideal_lattice)
    .def_readwrite("assume_ideal_structure",&mapping::MappingInput::assume_ideal_structure)
    .def_readwrite("options",&mapping::MappingInput::options);
    }

    // clang-format off
    m.def("structure_score", &mapping::structure_score);
    m.def("map_structure", &mapping::map_structure);
    // clang-format on
}
} // namespace wrappy
