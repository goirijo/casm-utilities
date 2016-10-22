#include "lib/completer/handlers.hpp"
#include <casm/completer/Handlers.hh>

namespace casmUtilities
{
    UtilityHandler::UtilityHandler(const std::string &init_utility_tag, const std::function<LaunchRuleList(po::options_description&)> &init_initializer):
        m_tag(init_utility_tag)
    {
        m_argument_rules=init_initializer(m_desc);
    }

    po::options_description &UtilityHandler::desc()
    {
        return m_desc;
    }

    const std::string &UtilityHandler::tag() const
    {
        return m_tag;
    }

    const LaunchRuleList &UtilityHandler::argument_rules() const
    {
        return m_argument_rules;
    }

    namespace utilityProgramOptions
    {
        void add_help_suboption(po::options_description &handler_desc)
        {
            handler_desc.add_options()
                ("help,h", "Print list of available options for this utility");

            return;
        }

        void add_output_suboption(po::options_description &handler_desc)
        {
            handler_desc.add_options()
                ("output,o", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path()), "Target output file");

            return;
        }
    }

}
