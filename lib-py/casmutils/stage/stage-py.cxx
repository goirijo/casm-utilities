#include "casmutils/stage.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace WrapPy
{
    std::string test()
    {
        return "This was a test. Did it work?";
    }
}

PYBIND11_MODULE(_stage, m)
{
    using namespace pybind11;

    m.doc() = "Temporary module for testing out work in progress.";

    m.def("test", WrapPy::test);
}
