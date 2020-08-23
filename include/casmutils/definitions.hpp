#ifndef CASM_UTILS_DEFINITIONS
#define CASM_UTILS_DEFINITIONS

#include <casm/misc/CASM_math.hh>
#include <casm/misc/CASM_Eigen_math.hh>
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include <filesystem>

namespace casmutils
{
namespace fs = std::filesystem;
/* using nlohmann::json; */
using CASM::almost_equal;
using CASM::iround;
using CASM::lround;
using CASM::round;
} // namespace casmutils

namespace utilities
{
namespace fs = std::filesystem;
} // namespace utilities

namespace simplicity
{
namespace fs = casmutils::fs;
} // namespace simplicity

#endif
