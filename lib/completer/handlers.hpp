#ifndef CASM_UTILS_HANDLERS
#define CASM_UTILS_HANDLERS

#include <casm/CASM_global_definitions.hh>
#include <casm/completer/Handlers.hh>
#include <tuple>
#include <utility>
#include <functional>
#include <map>

/**
 * An extension of CASM::OptionHandlerBase, repurposed
 * for bash completion functionality of casm-utilities
 * instead of casm.
 */

namespace casmUtilities
{

    /**
     * Container class to store all the different kinds
     * of possible properties a UtilityHandler might
     * have to deal with. The number of types shouldn't be
     * too outrageous, so the list of members won't grow too
     * much. It basically works like std::map, but for multiple
     * types, defined within the class.
     *
     * Values are accessed through strings, since each property
     * is meant to be associated with a program option. If there
     * is not property associated with the given string, a new
     * one is added. Addition of new properties can be suppressed
     * via the ::lock() routine, which will ensure that
     * an assertion is used before returning the property.
     *
     * Properties are always returned by reference.
     */

    class UtilityHandlerProperties
    {
        public:

            typedef int Index;

            ///When locked, new properties cannot be added (works via an assert statement!)
            void lock() const;

            ///When unlocked, new properties can be added (works via an assert statement!)
            void unlock() const;

            ///Returns a reference to the desired property, which is held in one of the several vectors of this class
            template<typename T>
                T &find(const std::string &key);


        private:

            ///Given a particular program option tag, get the index into the appropriate vector that holds the associated property
            std::map<std::string, Index> m_property_index_map;

            ///If true, an assertion will check to make sure the string identifier was 
            mutable bool m_locked;

            ///List of int values
            std::vector<int> m_int_props;

            ///List of string values
            std::vector<std::string> m_string_props;

            ///Primary templated method for retrieving types
            template<typename T>
                T &_retreive(const std::string &key);
    };

    /**
     * An alternative to CASM::OptionHandlerBase. Does not require
     * any derived classes. Instead, values are passed in the form
     * of std::tuple and std::function that specifies how to construct
     * the options description.
     */

    /*
    template<typename... Params>
    class UtilityHandler
    {
        public:

            ///Define the name of the command during construction
            UtilityHandler(const std::string &init_utility_tag,
                    const std::function<void(std::tuple<Params...>&, CASM::po::option_description&)> &initializer);

        private:

            ///The desired name for the casm-utilities option
            const std::string m_tag;

            ///Boost program options, which will hold all the available --suboptions for the utility
            CASM::po::option_description m_desc;

            ///Holds all the parameters whose values are filled through command line parsing
            std::tuple<Params...> m_parameters;
    };
    */

    class SplashOption : public CASM::Completer::OptionHandlerBase
    {
        public:

            SplashOption();

            int number() const;

        private:

            void initialize() override;

            int m_number;
    };

    class StructureOption : public CASM::Completer::OptionHandlerBase
    {
        public:

            using CASM::Completer::OptionHandlerBase::output_path;

            StructureOption();

        private:

            void initialize() override;
    };

}

#endif

