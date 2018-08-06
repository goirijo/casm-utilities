#include <iostream>

#include "casmutils/exceptions.hpp"
#include "casmutils/handlers.hpp"
#include <casm/completer/Handlers.hh>

namespace Utilities
{
Handler::Handler(int argc, char* argv[], const std::function<void(po::options_description&)>& init_initializer)
    : m_argc(argc)
{
    init_initializer(m_desc);
    po::store(po::parse_command_line(argc, argv, m_desc), m_vm);

    /* if (!m_utility.argument_rules().parse(m_vm)) */
    /* { */
    /*     throw UserInputMangle("A forbidden combination of command line arguments was parsed."); */
    /* } */
}

bool Handler::count(const std::string& countable) const { return m_vm.count(countable); }

void Handler::notify() { po::notify(m_vm); }

const po::variables_map& Handler::vm() const { return m_vm; }

int Handler::argc() const { return m_argc; }

const po::options_description& Handler::desc() const { return m_desc; }

namespace UtilityProgramOptions
{
void add_help_suboption(po::options_description& handler_desc)
{
    handler_desc.add_options()("help,h", "Print list of available options for this utility");

    return;
}

void add_desc_suboption(po::options_description& handler_desc)
{
    handler_desc.add_options()("desc", "Extended usage description");

    return;
}
void add_output_suboption(po::options_description& handler_desc)
{
    handler_desc.add_options()("output,o", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path()),
                               "Target output file");

    return;
}
} // namespace UtilityProgramOptions
} // namespace Utilities
