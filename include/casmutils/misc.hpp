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
 * of casm-utilities (e.g. Simplicity).
 * Usage of anything within this namespace should not leak
 * outside of any implementation, these are not utility
 * library functions, they're just here for convenience.
 */

namespace Extend
{
/// Make a PrimClex but shut up about it
CASM::PrimClex quiet_primclex(CASM::Structure& prim);
} // namespace extend

#endif
