#ifndef CASM_UTILS_LAUNCH
#define CASM_UTILS_LAUNCH

#include "lib/completer/handlers.hpp"
#include "lib/definitions.hpp"
#include "lib/launch/rules.hpp"

namespace casmUtilities
{
/**
 * Container class for parsing command line
 * arguments and filling up all the boost
 * variables needed
 */

class Launcher
{
public:
    Launcher(int argc,
             char *argv[],
             const std::string &init_utility_tag,
             const std::function<LaunchRuleList(po::options_description &)> &init_initializer);

    /// Check to see if the given string was parsed into the variables map
    bool count(const std::string &countable) const;

    UtilityHandler &utility();

    void notify();

    const po::variables_map &vm() const;

    template <typename T>
    T fetch(const std::string &prog_option) const;

    int argc() const;

private:
    /// Keep info about raw command line input (number of arguments)
    int m_argc;

    /// Derived from CASM::OptionHandlerBase
    UtilityHandler m_utility;

    /// Deals with the given command line options
    po::variables_map m_vm;
};

//************************************************************************************//

template <typename T>
T Launcher::fetch(const std::string &prog_option) const
{
    return vm()[prog_option].as<T>();
}
}

#endif
