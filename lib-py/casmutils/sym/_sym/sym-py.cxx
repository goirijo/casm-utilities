#include <casmutils/sym/cartesian.hpp>
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
PYBIND11_MODULE(_sym, m)
{
    using namespace pybind11;

    m.doc() = "Python bindings for classes and functions related symmetry representations and groups\
               , such as crystal structures or clusters.";

    {
        //TODO: Have a copy of this class with only def_read and then make this one MutableCartOp?
        //or should we just stick to that distinction for classes that have non-const functions?
        class_<sym::CartOp>(m, "CartOp")
            .def(init<const Eigen::Matrix3d&, const Eigen::Vector3d&, bool>())
            .def_static("identity",&sym::CartOp::identity)
            .def_static("time_reversal",&sym::CartOp::time_reversal)
            .def_static("translation_operation",&sym::CartOp::translation_operation)
            .def_static("point_operation",&sym::CartOp::point_operation)
            .def_readwrite("matrix",&sym::CartOp::matrix)
            .def_readwrite("translation",&sym::CartOp::translation)
            .def_readwrite("is_time_reversal_active",&sym::CartOp::is_time_reversal_active)
            .def("__mul__",(sym::CartOp(*)(const sym::CartOp&, const sym::CartOp&))sym::operator*);
    }
}
} // namespace wrappy
