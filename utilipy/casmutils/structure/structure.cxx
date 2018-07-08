//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/operators.hpp>
#include <boost/filesystem.hpp>

#include <casm/CASM_global_definitions.hh>
#include "casmutils/structure.hpp"
#include <casm/crystallography/Structure.hh>
#include <string>

char const* greet()
{
   return "hello, world";
}

std::string smart_greet()
{
    std::string greeting="Can you do this though?";
    return std::string();
}

namespace Rewrap
{
    class Structure : public CASM::Structure
    {
        public:

            Structure()=delete;

            Structure(CASM::Structure init_struc):CASM::Structure(init_struc){}

            bool is_primitive() const
            {
                return CASM::Structure::is_primitive();
            }

            Structure primitive() const
            {
                return Simplicity::make_primitive(*this);
            }

        private:


    };

    Structure from_poscar(const std::string& filename)
    {
        return Structure(CASM::Structure(filename));
    }

    void to_poscar(const Structure& writeable, const std::string& filename)
    {
        Simplicity::write_poscar(writeable, filename);
        return;
    }

    std::ostream& operator<<(std::ostream& output, const Structure& outputtable)
    {
        Simplicity::print_poscar(outputtable, output);
        return output;
    }

}

BOOST_PYTHON_MODULE(_structure)
{
    using namespace boost::python;
    using namespace Rewrap;

    class_<Structure>("Structure",no_init)
        .def("is_primitive", &Structure::is_primitive)
        .def("primitive", &Structure::primitive)
        .def(str(self))
        ;

    def("from_poscar",from_poscar);
    def("to_poscar",to_poscar);
}


