/* #include <iostream> */

/* #include "casm/completer/Handlers.hh" */
/* #include "casmutils/exceptions.hpp" */
/* #include "casmutils/handlers.hpp" */

/* namespace utilities */
/* { */
/* Handler::Handler(int argc, char* argv[], const std::function<void(po::options_description&)>& init_initializer) */
/*     : m_argc(argc) */
/* { */
/*     init_initializer(m_desc); */
/*     po::store(po::parse_command_line(argc, argv, m_desc), m_vm); */

/*     /1* if (!m_utility.argument_rules().parse(m_vm)) *1/ */
/*     /1* { *1/ */
/*     /1*     throw UserInputMangle("A forbidden combination of command line arguments was parsed."); *1/ */
/*     /1* } *1/ */
/* } */

/* bool Handler::count(const std::string& countable) const { return m_vm.count(countable); } */

/* void Handler::notify() { po::notify(m_vm); } */

/* const po::variables_map& Handler::vm() const { return m_vm; } */

/* int Handler::argc() const { return m_argc; } */

/* const po::options_description& Handler::desc() const { return m_desc; } */

/* namespace utilities */
/* { */
/* void add_help_suboption(po::options_description& handler_desc) */
/* { */
/*     handler_desc.add_options()("help,h", "Print list of available options for this utility"); */

/*     return; */
/* } */

/* void add_desc_suboption(po::options_description& handler_desc) */
/* { */
/*     handler_desc.add_options()("desc", "extended usage description"); */

/*     return; */
/* } */
/* void add_output_suboption(po::options_description& handler_desc) */
/* { */
/*     handler_desc.add_options()("output,o", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path()), */
/*                                "Target output file"); */

/*     return; */
/* } */
/* } // namespace utilities */
/* } // namespace utilities */
