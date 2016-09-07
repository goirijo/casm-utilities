#include <casm/CASM_global_definitions.hh>

#include "lib/completer/handlers.hpp"

namespace casmUtilities
{

    template<>
        int &UtilityHandlerProperties::_retreive(const std::string &key)
        {
            return m_int_props[m_property_index_map[key]];
        }

    template<>
        std::string &UtilityHandlerProperties::_retreive(const std::string &key)
        {
            return m_string_props[m_property_index_map[key]];
        }

    template<typename T>
        T &UtilityHandlerProperties::find(const std::string &key)
        {
            if(m_locked)
            {
                assert(m_property_index_map.find(key)!=m_property_index_map.end());
            }

            return _retreive<T>(key);
        }

    /*
    UtilityHandler::UtilityHandler(const std::string &init_utility_tag,
                    const std::function<void(std::tuple<Params...>&, CASM::po::option_description&)> &initializer)
    {
    }
    */

    SplashOption::SplashOption() : OptionHandlerBase("splash", "A simple command that does nothing useful"){}

    void SplashOption::initialize()
    {
        add_help_suboption();

        m_desc.add_options()
            ("print,p", "Print the CASM logo to the screen")
            ("dumb-print,d", "Print the CASM logo as a word")
            ("squelch,s", "Don't print anything")
            ("test,t", CASM::po::value<int>(), "Test without storage")
            ("number,n", CASM::po::value<int>(&m_number)->default_value(1), "How many times to print the logo");


        return;
    }

    int SplashOption::number() const
    {
        return m_number;
    }

    StructureOption::StructureOption() : OptionHandlerBase("structure", "Manipulate your basis and unit cells"){}

    void StructureOption::initialize()
    {
        add_help_suboption();

        return;
    }

}
