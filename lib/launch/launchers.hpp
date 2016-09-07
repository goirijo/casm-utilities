#ifndef CASM_UTILS_LAUNCH
#define CASM_UTILS_LAUNCH

#include <casm/CASM_global_definitions.hh>
#include "lib/launch/rules.hpp"

namespace casmUtilities
{
    /**
     * Container class for parsing command line
     * arguments and filling up all the boost
     * variables needed
     */

    template <typename T>
    class Launcher
    {
        public:

            Launcher(int argc, char *argv[]);

            bool count(const std::string &countable) const;

            T &option();

            const CASM::po::variables_map &vm() const;

        private:

            ///Derived from CASM::OptionHandlerBase
            T m_option;

            ///Deals with the given command line options
            CASM::po::variables_map m_vm;

            ///Specify how the command line options should interact with each other
            LaunchRuleList m_argument_rules;
    };

    //************************************************************************************//

    template<typename T>
    Launcher<T>::Launcher(int argc, char *argv[])
    {
        CASM::po::store(CASM::po::parse_command_line(argc, argv, m_option.desc()),m_vm);
        CASM::po::notify(m_vm);
    }

    template<typename T>
    bool Launcher<T>::count(const std::string &countable) const
    {
        return m_vm.count(countable);
    }

    template<typename T>
    T &Launcher<T>::option()
    {
        return m_option;
    }

    template<typename T>
    const CASM::po::variables_map &Launcher<T>::vm() const
    {
        return m_vm;
    }


}

#endif
