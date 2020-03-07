#ifndef SITE_PY_HH
#define SITE_PY_HH

#include <string>

namespace casmutils::xtal
{
class Site;
}

namespace wrappy
{
using namespace casmutils;
namespace Site
{
std::string __str__(const xtal::Site& printable);
} // namespace Site
} // namespace wrappy

#endif
