#ifndef CASM_UTILS_DEFINITIONS
#define CASM_UTILS_DEFINITIONS

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"

namespace Rewrap
{
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using CASM::COORD_TYPE;
using CASM::CART;
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
