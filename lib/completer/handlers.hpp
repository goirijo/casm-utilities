#ifndef CASM_UTILS_HANDLERS
#define CASM_UTILS_HANDLERS

#include <casm/CASM_global_definitions.hh>
#include <casm/completer/Handlers.hh>
#include <utility>
#include <functional>

#include "lib/definitions.hpp"
#include "lib/launch/rules.hpp"

/**
 * An extension of CASM::OptionHandlerBase, repurposed
 * for bash completion functionality of casm-utilities
 * instead of casm.
 */

namespace casmUtilities
{

    /**
     * An alternative to CASM::OptionHandlerBase. Does not require
     * any derived classes. Instead, an std::function is passed that
     * specifies how to construct the program options. The only restriction
     * for this is that none of the values that get parsed from the
     * command line are allowed to be stored in predefined variables.
     * Don't use something like po::value<int>(&my_int), just do
     * po::value<int>. Internally, vm.as<int>() will be used to
     * access the value.
     */

    class UtilityHandler
    {
        public:

            ///Define the name of the command during construction
            UtilityHandler(const std::string &init_utility_tag, const std::function<LaunchRuleList(po::options_description&)> &init_initializer);
            
            ///Get the program options, filled with initialized values
            po::options_description &desc();

            ///The desired name for the utility
            const std::string &tag() const;

        private:

            ///The desired name for the casm-utilities option
            const std::string m_tag;

            ///Boost program options, which will hold all the available --suboptions for the utility
            po::options_description m_desc;

            ///Specify how the command line options should interact with each other
            LaunchRuleList m_argument_rules;

    };

}

#endif

