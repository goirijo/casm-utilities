#include <boost/filesystem.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casmutils/structure.hpp"
#include <casm/crystallography/Structure.hh>
#include <string>

//******************************************************************************************************//
//******************************************************************************************************//

namespace WrapPy
{

PYBIND11_MODULE(_frankenstein, m)
{
    using namespace pybind11;
    using namespace WrapPy;

    m.doc() = "Structure manipulation for slicing, shifting, and stacking.";

    m.def("from_poscar", from_poscar);
    m.def("to_poscar", to_poscar);
    /* m.def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&)) Simplicity::make_niggli); */
}
}
