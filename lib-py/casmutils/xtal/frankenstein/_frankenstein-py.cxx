#include <casm/external/Eigen/Dense>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casmutils/xtal/frankenstein.hpp"
#include <string>

//******************************************************************************************************//
//******************************************************************************************************//

namespace wrappy
{

rewrap::Structure shift_coords_by(rewrap::Structure* struc, const Eigen::Vector3d& shift_val)
{
    auto mutable_struc = *struc;
    frankenstein::shift_coords_by(&mutable_struc, shift_val);
    return mutable_struc;
}

PYBIND11_MODULE(_frankenstein, m)
{
    using namespace pybind11;
    using namespace wrappy;

    m.doc() = "Structure manipulation for slicing, shifting, and stacking.";

    m.def("slice", frankenstein::slice);
    m.def("multi_slice", frankenstein::multi_slice);
    m.def("uniformly_slice", frankenstein::uniformly_slice);
    m.def("stack", frankenstein::stack);
    m.def("vacuum_pack", frankenstein::vacuum_pack);
    m.def("inflate", frankenstein::inflate);
    m.def("shift_coords_by", shift_coords_by);
}
} // namespace wrappy
