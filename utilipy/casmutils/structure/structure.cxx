//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/filesystem.hpp>
/* #include <boost/python/class.hpp> */
/* #include <boost/python/def.hpp> */
/* #include <boost/python/list.hpp> */
/* #include <boost/python/module.hpp> */
/* #include <boost/python/operators.hpp> */
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "casmutils/structure.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <string>

#include "Eigen/Core"

//******************************************************************************************************//
//******************************************************************************************************//

namespace Rewrap
{

class Structure : public CASM::Structure
{
public:
    Structure() = delete;

    Structure(CASM::Structure init_struc);

    bool is_primitive() const;

    Structure primitive() const;

    std::string to_string() const;

private:
};

Structure from_poscar(const std::string& filename);
void to_poscar(const Structure& writeable, const std::string& filename);

//******************************************************************************************************//

Structure::Structure(CASM::Structure init_struc) : CASM::Structure(init_struc) {}

bool Structure::is_primitive() const { return CASM::Structure::is_primitive(); }

Structure Structure::primitive() const { return Simplicity::make_primitive(*this); }

std::string Structure::to_string() const
{
    std::ostringstream sstream;
    Simplicity::print_poscar(*this, sstream);
    return sstream.str();
}

Structure make_niggli(const Structure& non_niggli)
{
    return Simplicity::make_niggli(non_niggli);
}

Structure from_poscar(const std::string& filename) { return Structure(CASM::Structure(filename)); }

void to_poscar(const Structure& writeable, const std::string& filename)
{
    Simplicity::write_poscar(writeable, filename);
    return;
}

void test_vec_int(const std::vector<int>& vec)
{
    for(const auto& v : vec)
    {
        std::cout<<v<<std::endl;
    }
}

void test_eigen_mat(const Eigen::Matrix3d& mat)
{
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            std::cout<<mat(i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    return;
}

/* std::ostream& operator<<(std::ostream& output, const Structure& outputtable) */
/* { */
/*     Simplicity::print_poscar(outputtable, output); */
/*     return output; */
/* } */

}

//******************************************************************************************************//

/* BOOST_PYTHON_MODULE(_structure) */
/* { */
/*     using namespace boost::python; */
/*     using namespace Rewrap; */

/*     class_<Structure>("Structure", no_init) */
/*         .def("is_primitive", &Structure::is_primitive) */
/*         .def("primitive", &Structure::primitive) */
/*         .def(str(self)); */

/*     /1* typedef std::vector<int> std_vector_int; *1/ */
/*     /1* class_<std_vector_int>("std_vector_int").def(vector_indexing_suite<std_vector_int>()); *1/ */

/*     def("from_poscar", from_poscar); */
/*     def("to_poscar", to_poscar); */
/*     /1* def("test_int_list", test_int_list); *1/ */
/* } */

PYBIND11_MODULE(_structure, m)
{
    using namespace pybind11;
    using namespace Rewrap;

    m.doc() = "Raw python bindings for a re-wrapped CASM::Structure class.";

    class_<Structure>(m,"Structure")
        .def("is_primitive", &Structure::is_primitive)
        .def("primitive", &Structure::primitive)
        .def("__repr__", &Structure::to_string)
        ;

    m.def("from_poscar", from_poscar);
    m.def("to_poscar", to_poscar);
    m.def("make_niggli", make_niggli);
    /* m.def("test_vec_int", test_vec_int); */
    /* m.def("test_eigen_mat", test_eigen_mat); */
}
