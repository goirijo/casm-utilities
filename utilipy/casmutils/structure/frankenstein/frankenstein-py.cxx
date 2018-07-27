#include <boost/filesystem.hpp>
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casmutils/frankenstein.hpp"
#include <casm/crystallography/Structure.hh>
#include <string>

//******************************************************************************************************//
//******************************************************************************************************//

namespace WrapPy
{

Rewrap::Structure shift_coords_by(Rewrap::Structure* struc, const Eigen::Vector3d& shift_val)
{
    auto mutable_struc = *struc;
    Frankenstein::shift_coords_by(&mutable_struc, shift_val);
    return mutable_struc;
}

PYBIND11_MODULE(_frankenstein, m)
{
    using namespace pybind11;
    using namespace WrapPy;

    m.doc() = "Structure manipulation for slicing, shifting, and stacking.";

    m.def("slice", Frankenstein::slice);
    m.def("multi_slice", Frankenstein::multi_slice);
    m.def("uniformly_slice", Frankenstein::uniformly_slice);
    m.def("stack", Frankenstein::stack);
    m.def("vacuum_pack", Frankenstein::vacuum_pack);
    m.def("inflate", Frankenstein::inflate);
    m.def("shift_coords_by", shift_coords_by);
}
} // namespace WrapPy
