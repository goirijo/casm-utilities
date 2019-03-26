#include <casm/clex/PrimClex.hh>

namespace Extend
{
CASM::PrimClex quiet_primclex(CASM::Structure& prim)
{
    CASM::Log log(std::cout, 0);
    CASM::Logging logging(log);
    CASM::PrimClex pclex(prim, logging);
    return pclex;
}
} // namespace extend
