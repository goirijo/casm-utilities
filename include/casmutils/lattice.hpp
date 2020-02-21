#ifndef UTILS_LATTICE_HH
#define UTILS_LATTICE_HH

#include "casmutils/definitions.hpp"
#include <casm/crystallography/Lattice.hh>

namespace CASM
{
}

namespace Rewrap
{
    class Lattice
    {
        public:

            Lattice() = delete;
            Lattice(const CASM::xtal::Lattice& init_lat);
            Lattice(const Eigen::Matrix3d& column_lat_mat);

            //TODO: Read up on Eigen::Matrix3d::ColXpr and decide if you prefer this. CASM does it this way.
            /// Return the ith vector of the lattice
            Eigen::Vector3d operator[](int i) const
            {
                return this->__get()[i];
            }

            /// Return the first vector of the lattice
            Eigen::Vector3d a() const
            {
                return this->operator[](0);
            }

            /// Return the second vector of the lattice
            Eigen::Vector3d b() const
            {
                return this->operator[](1);
            }

            /// Return the third vector of the lattice
            Eigen::Vector3d c() const
            {
                return this->operator[](2);
            }

            double volume() const
            {
                return this->__get().volume();
            }
            Eigen::Matrix3d column_vector_matrix() const {
                return this->__get().lat_column_mat();
            }
            ///Return *this as a CASM::Lattice
            const CASM::xtal::Lattice& __get() const {return this->casm_lattice;}

        private:

            CASM::xtal::Lattice casm_lattice;
    };
}

#endif
