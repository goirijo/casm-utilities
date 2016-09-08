#include <iostream>
#include <iterator>

#include "lib/launch/launchers.hpp"

namespace casmUtilities
{
    Launcher::Launcher(int argc, 
            char *argv[],
            const std::string &init_utility_tag,
            const std::function<LaunchRuleList(po::options_description&)> &init_initializer):
        m_utility(init_utility_tag, init_initializer)
    {
        po::store(po::parse_command_line(argc, argv, m_utility.desc()), m_vm);
        po::notify(m_vm);

        if(!m_utility.argument_rules().parse(m_vm))
        {
            throw UserInputMangle("A forbidden combination of command line arguments was parsed.");
        }
    }

    bool Launcher::count(const std::string &countable) const
    {
        return m_vm.count(countable);
    }

    UtilityHandler &Launcher::utility()
    {
        return m_utility;
    }

    const po::variables_map &Launcher::vm() const
    {
        return m_vm;
    }

}
