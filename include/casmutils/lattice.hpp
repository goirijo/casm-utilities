#ifndef UTILS_LATTICE_HH
#define UTILS_LATTICE_HH

#include "casmutils/definitions.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Lattice.hh>

namespace CASM
{
}

namespace Rewrap
{
    class Lattice : public CASM::Lattice
    {
        public:

            Lattice() = delete;
            Lattice(const CASM::Lattice& init_lat);
            Lattice(const Eigen::Matrix3d& column_lat_mat);

            ///Return *this as a CASM::Lattice
            const CASM::Lattice& __get() const {return *this;};

        private:
    };
}

#endif
