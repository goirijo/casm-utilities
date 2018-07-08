//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
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
    return std::string("Can you do this though?");
}

namespace Rewrap
{
    class Structure
    {
        public:

            Structure(const std::string& poscar_name): m_core(CASM::fs::path(poscar_name))
            {
            }

            const CASM::Structure& get() const
            {
                return m_core;
            }

        private:

            CASM::Structure m_core;
    };

    void write_poscar(const Structure& printable, const std::string& filename)
    {
        return Simplicity::write_poscar(printable.get(), CASM::fs::path(filename));
    }
}

BOOST_PYTHON_MODULE(hello_ext)
{
    using namespace boost::python;
    def("greet", greet);
    def("smart_greet", smart_greet);

    using namespace boost::python;
    using namespace Rewrap;

    class_<Structure>("Structure",init<std::string>());

    def("write_poscar", write_poscar);
}


