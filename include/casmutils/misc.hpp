#ifndef MISC_HH
#define MISC_HH

namespace CASM
{
class PrimClex;
class Structure;
} // namespace CASM

/**
 * This namespace is reserved for convenience functions
 * that reduce boilerplate code within library functions
 * of casm-utilities (e.g. Simplicity)
 */

namespace extend
{
/// Make a PrimClex but shut up about it
CASM::PrimClex quiet_primclex(CASM::Structure& prim);
} // namespace extend

#endif
