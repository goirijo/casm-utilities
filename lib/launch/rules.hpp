#ifndef CASM_UTILS_RULES
#define CASM_UTILS_RULES

#include <vector>
#include <string>
#include <utility>
#include <casm/CASM_global_definitions.hh>

namespace casmUtilities
{
    /**
     * Defines whether the amount of --suboptions given to a particular
     * utility was underspecified (requires --more input), overspecified
     * (--conflicting --options were given), correctly specified, or
     * unspecified (no rule given).
     */

    enum class RuleStatus {UNDERSPECIFIED, OVERSPECIFIED, SPECIFIED, UNSPECIFIED};

    /**
     * After you've parsed the command line and constructed a po::variables_map,
     * you need to ensure that conflicting arguments haven't been given. This
     * class allows you to specify which --suboptions conflict with others, and
     * which ones must be simultaneously specified. After building up
     * an instance you can then check to see whether the variables map you
     * have satisfies all the requirements.
     */

    class LaunchRule
    {
        public:

            /// \brief specify single --suboption that cannot go with another
            void add_exclusion(std::string &target, std::string &exclude);

            /// \brief specify single --suboption that must go with another
            void add_inclusion(std::string &target, std::string &include);

            /// \brief given a collection of --suboptions, require that if one of them is given, all of them must be given
            void add_requirement(const std::vector<std::string> &codependent_group);

            /// \brief given a --suboption, require non of the other --suboptions in the other group are given
            void add_constraint(const std::string &independent, const std::vector<std::string> &exclusion_group);

            /// \brief Using the loaded rules, determine whether the given variables map is good to go
            RuleStatus parse(const CASM::po::variables_map &ref_vm);

        private:

            /// \brief List of every --suboption that excludes others
            std::vector<std::pair<std::string, std::string> > m_exclusion_table;
            
            /// \brief List of every --suboption that requires others
            std::vector<std::pair<std::string, std::string> > m_inclusion_table;

            /// \brief List of every rule broken since the last ::parse
            std::vector<std::pair<std::string, std::string> > m_broken_rules;
    };
}


#endif
