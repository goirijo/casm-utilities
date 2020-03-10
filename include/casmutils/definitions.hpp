#ifndef CASM_UTILS_DEFINITIONS
#define CASM_UTILS_DEFINITIONS

#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include <boost/program_options.hpp>
#include <filesystem>

namespace casmutils
{
namespace fs = std::filesystem;
namespace po = boost::program_options;
} // namespace casmutils

namespace utilities
{
namespace fs = std::filesystem;
namespace po = boost::program_options;
} // namespace utilities

namespace simplicity
{
namespace fs = casmutils::fs;
namespace po = casmutils::po;
} // namespace simplicity

#endif
