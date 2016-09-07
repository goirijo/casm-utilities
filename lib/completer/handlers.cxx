#include <casm/CASM_global_definitions.hh>

#include "lib/completer/handlers.hpp"

namespace casmUtilities
{
    UtilityHandler::UtilityHandler(const std::string &init_utility_tag, const std::function<LaunchRuleList(po::options_description&)> &init_initializer)
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

}
