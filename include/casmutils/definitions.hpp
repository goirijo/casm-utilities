#ifndef CASM_UTILS_DEFINITIONS
#define CASM_UTILS_DEFINITIONS

#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace Rewrap
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using CASM::CART;
using CASM::COORD_TYPE;
using CASM::FRAC;
} // namespace Rewrap

namespace Simplicity
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;
} // namespace Simplicity

namespace Utilities
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;
} // namespace Utilities

#endif
