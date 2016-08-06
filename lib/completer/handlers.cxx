#include <casm/CASM_global_definitions.hh>

#include "lib/completer/handlers.hpp"

namespace casmUtilitiesCompletion
{
    SplashOption::SplashOption() : OptionHandlerBase("splash", "A simple command that does nothing useful"){}

    void SplashOption::initialize()
    {
        add_help_suboption();

        m_desc.add_options()
            ("print,p", "Print the CASM logo to the screen")
            ("dumb-print,d", "Print the CASM logo as a word")
            ("squelch,s", "Don't print anything")
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
