#ifndef CASM_UTILS_LAUNCH
#define CASM_UTILS_LAUNCH

#include "casmutils/definitions.hpp"
/* #include "casmutils/handlers.hpp" */
/* #include "casmutils/rules.hpp" */

namespace Utilities
{
    /**
     * Container class for parsing command line
     * arguments and filling up all the boost
     * variables needed
     */

    class Handler
    {
        public:

            Handler(int argc, 
                    char *argv[],
                    const std::function<void(po::options_description&)> &init_initializer);

            ///Check to see if the given string was parsed into the variables map
            bool count(const std::string &countable) const;

            /* UtilityHandler &utility(); */

            void notify();

            const po::variables_map &vm() const;

            template<typename T>
                T fetch(const std::string &prog_option) const;

             int argc() const;

             const po::options_description& desc() const;

        private:

            /// Keep info about raw command line input (number of arguments)
            int m_argc;

            ///Boost program options, which will hold all the available --suboptions for the utility
            po::options_description m_desc;

            ///Deals with the given command line options
            po::variables_map m_vm;

    };

    //************************************************************************************//
    
    template<typename T>
    T Handler::fetch(const std::string &prog_option) const
    {
        return vm()[prog_option].as<T>();
    }
    
    //************************************************************************************//

    /**
     * A bunch of commonly used program options to reuse
     * across different UtilityHandler objects
     */

    namespace UtilityProgramOptions
    {
        /// \brief Lame help message with each available option
        void add_help_suboption(po::options_description &handler_desc);

        /// \brief Lame help message with each available option
        void add_desc_suboption(po::options_description &handler_desc);

        /// \brief Store a fs::path to which whatever you're going to output should be saved
        void add_output_suboption(po::options_description &handler_desc);
    }


}

#endif
