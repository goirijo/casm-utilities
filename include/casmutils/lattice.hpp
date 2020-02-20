#ifndef UTILS_LATTICE_HH
#define UTILS_LATTICE_HH

#include "casmutils/definitions.hpp"
#include <casm/crystallography/Lattice.hh>

namespace CASM
{
}

namespace Rewrap
{
    class Lattice : public CASM::xtal::Lattice
    {
        public:

            Lattice() = delete;
            Lattice(const CASM::xtal::Lattice& init_lat);
            Lattice(const Eigen::Matrix3d& column_lat_mat);

            ///Return *this as a CASM::Lattice
            const CASM::xtal::Lattice& __get() const {return *this;};

        private:
    };
}

#endif
