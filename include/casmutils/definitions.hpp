#ifndef CASM_UTILS_DEFINITIONS
#define CASM_UTILS_DEFINITIONS

#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include <filesystem>
#include <boost/program_options.hpp>

namespace CASMUtils
{
namespace fs = std::filesystem;
namespace po = boost::program_options;
} // namespace Utilities

namespace Utilities
{
namespace fs = std::filesystem;
namespace po = boost::program_options;
} // namespace Utilities

namespace rewrap
{
namespace fs=CASMUtils::fs;
namespace po=CASMUtils::po;
using CASM::CART;
using CASM::COORD_TYPE;
using CASM::FRAC;
} // namespace rewrap

namespace Simplicity
{
namespace fs=CASMUtils::fs;
namespace po=CASMUtils::po;
} // namespace Simplicity

#endif
