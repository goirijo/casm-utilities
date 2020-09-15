#include <casmutils/sym/cartesian.hpp>
#include <fstream>
#include <string>

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
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
PYBIND11_MODULE(_sym, m)
{
    using namespace pybind11;

    m.doc() = "Python bindings for classes and functions related symmetry representations and groups\
               , such as crystal structures or clusters.";

    {
        class_<sym::CartOp>(m, "CartOp")
            .def(init<const Eigen::Matrix3d&, const Eigen::Vector3d&, bool>())
            .def(init<const sym::CartOp&>())
            .def_static("identity", &sym::CartOp::identity)
            .def_static("time_reversal", &sym::CartOp::time_reversal)
            .def_static("translation_operation", &sym::CartOp::translation_operation)
            .def_static("point_operation", &sym::CartOp::point_operation)
            .def_readwrite("matrix", &sym::CartOp::matrix)
            .def_readwrite("translation", &sym::CartOp::translation)
            .def_readwrite("is_time_reversal_active", &sym::CartOp::is_time_reversal_active)
            .def(pybind11::self * pybind11::self);
    }
}
} // namespace wrappy
