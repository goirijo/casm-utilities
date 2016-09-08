#ifndef CASM_UTILS_RULES
#define CASM_UTILS_RULES

#include <vector>
#include <string>
#include <utility>

#include "lib/definitions.hpp"

namespace casmUtilities
{
    /**
     * Defines whether the amount of --suboptions given to a particular
     * utility was underspecified (requires --more input), overspecified
     * (--conflicting --options were given), correctly specified, or
     * unspecified (no rule given).

    enum class RuleStatus {UNDERSPECIFIED, OVERSPECIFIED, SPECIFIED, UNSPECIFIED};
     */

    /**
     * Defines what kind of a rule the --suboption pair defines:
     *  -If the first string is given the second must ALSO be given
     *  -If the first string is given the second must NOT be given
     *  -NO strings should be given
     *  -At least ONE string must be given
     */

    enum class RuleType {INCLUSIVE, EXCLUSIVE, SILENT, REQUIRED};

    /**
     * After you've parsed the command line and constructed a po::variables_map,
     * you need to ensure that conflicting arguments haven't been given. This
     * class allows you to specify which --suboptions conflict with others, and
     * which ones must be simultaneously specified.
     */

    class LaunchRule
    {
    public:

        /// \brief Explicitly declare everything needed
        LaunchRule(const std::string &target, const std::string &requirement, RuleType specification);

        /// \brief Only allowed for SILENT and REQUIRED
        LaunchRule(RuleType restricted_specification);

        /// \brief access what kind of rule *this is
        RuleType type() const;

        /// \brief check to see if the given variables map breaks the rule (false->invalid, true->valid)
        bool parse(const po::variables_map &ref_vm) const;

        /// \brief check the specific strings of the rule
        const std::pair<std::string, std::string> &requirement_pair() const;
        

    private:

        /// \brief Default constructor is private
        LaunchRule();

        /// \brief a pair of --suboptions, that are either required or exclude each other
        std::pair<std::string, std::string> m_requirement_pair;

        /// \brief defines the basic type of rule
        RuleType m_type;

        /// \brief parse variables map when *this is INCLUSIVE
        bool _inclusive_parse(const po::variables_map &ref_vm) const;

        /// \brief parse variables map when *this is EXCLUSIVE
        bool _exclusive_parse(const po::variables_map &ref_vm) const;

        /// \brief parse variables map when *this is SILENT
        bool _silent_parse(const po::variables_map &ref_vm) const;

        /// \brief parse variables map when *this is REQUIRED
        bool _required_parse(const po::variables_map &ref_vm) const;

    };

    /**
     * A collection of LaunchRules that can be parsed. After building up
     * an instance you can then check to see whether the variables map you
     * have satisfies all the requirements.
     */

    class LaunchRuleList
    {
    public:

        /// \brief specify single --suboption that cannot go with another
        void add_exclusion(const std::string &target, const std::string &exclude);

        /// \brief specify single --suboption that must go with another
        void add_inclusion(const std::string &target, const std::string &include);

        /// \brief specify single --suboption that must go with *only* one of the given other --suboptions
        void add_any_inclusion(const std::string &target, const std::vector<std::string> &includes);

        /// \brief specify single --suboption that can't go with *any* of the given other --suboptions
        void add_any_exclusion(const std::string &target, const std::vector<std::string> &exclude);



        /// \brief given a collection of --suboptions, require that if one of them is given, all of them must be given
        void add_requirement(const std::vector<std::string> &codependent_group);

        /// \brief given a --suboption, require none of the other --suboptions in the other group are given
        void add_constraint(const std::string &independent, const std::vector<std::string> &exclusion_group);

        /// \brief Using the loaded rules, determine whether the given variables map is good to go
        bool parse(const po::variables_map &ref_vm) const;

    private:

        /// \brief List of every rule that must be obeyed
        std::vector<LaunchRule> m_strict_rule_table;

        /// \brief List of groups of rules of which only one rule must be obeyed
        std::vector<std::vector<LaunchRule> > m_any_rule_table;

        /// \brief List of every rule broken since the last ::parse related to the strict set
        mutable std::vector<int> m_strict_broken_rules;

        /// \brief List of every rule broken since the last ::parse related to the any set
        mutable std::vector<int> m_any_broken_rules;
    };

    class UserInputMangle: public std::runtime_error
    {
        public:
            UserInputMangle(const std::string &init_message):std::runtime_error(init_message) {}
    };
}


#endif
