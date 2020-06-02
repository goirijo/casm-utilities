#ifndef UTILS_SYM_CART_HH
#define UTILS_SYM_CART_HH

#include <casm/crystallography/SymType.hh>
#include <variant>
#include <vector>

namespace casmutils
{
namespace sym
{
typedef CASM::xtal::SymOp CartOp;
typedef std::vector<int> PermRep;
using CASM::xtal::operator*;
} // namespace sym
} // namespace casmutils

#endif
