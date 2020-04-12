#include "casmutils/xtal/structure.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casm/crystallography/SymTools.hh>
#include <casm/crystallography/BasicStructureTools.hh>

namespace casmutils
{
    namespace xtal
    {
        std::vector<sym::CartOp> make_point_group(const Lattice& lat, double tol)
        {
            return CASM::xtal::make_point_group(lat.__get(), tol);
        }

        std::vector<sym::CartOp> make_factor_group(const Structure& struc, double tol)
        {
            return CASM::xtal::make_factor_group(struc.__get<CASM::xtal::BasicStructure>(),tol);
        }
    }
}
