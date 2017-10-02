#include<iostream>

#include "lib/launch/rules.hpp"

namespace casmUtilities
{

    LaunchRule::LaunchRule(const std::string &target, const std::string &requirement, RuleType specification):
        m_requirement_pair(std::make_pair(target, requirement)),
        m_type(specification)
    {}


    LaunchRule::LaunchRule(RuleType restricted_specification):
        m_type(restricted_specification)
    {
        assert(restricted_specification==RuleType::SILENT || restricted_specification==RuleType::REQUIRED);
    }

    RuleType LaunchRule::type() const
    {
        return m_type;
    }
    
    /**
     * Depending on the type of *this, the variables map will be parsed accordingly.
     * False is returned if the rule specification is broken, and true is returned
     * if the rule is met.
     */

    bool LaunchRule::parse(const po::variables_map &ref_vm) const
    {
        bool good_input=false;

        switch(m_type)
        {
            case RuleType::INCLUSIVE:
                good_input=_inclusive_parse(ref_vm);
                break;

            case RuleType::EXCLUSIVE:
                good_input=_exclusive_parse(ref_vm);
                break;

            default:
                std::cerr<<"UNIMPLEMENTED PARSE VALUE in rules.cxx"<<std::endl;
                assert(false);
        }

        return good_input;
    }

    const std::pair<std::string, std::string> &LaunchRule::requirement_pair() const
    {
        return m_requirement_pair;
    }

    /**
     * For the INCLUSIVE rule, if the first --suboption is given, then the
     * second suboption must *also* be given. If the first --suboption is not given,
     * the default parse value is returned.
     */

    bool LaunchRule::_inclusive_parse(const po::variables_map &ref_vm) const
    {
        assert(m_type==RuleType::INCLUSIVE);
        bool good_input;

        if(ref_vm.count(m_requirement_pair.first))
        {
            good_input=(ref_vm.count(m_requirement_pair.first) && ref_vm.count(m_requirement_pair.second));
        }

        else
        {
            good_input=true;
        }

        return good_input;
    }

    /**
     * For the EXCLUSIVE rule, if the first --suboption is given, then the
     * second suboption must *not* be given. If the first --suboption is not given,
     * the default parse value is returned.
     */

    bool LaunchRule::_exclusive_parse(const po::variables_map &ref_vm) const
    {
        assert(m_type==RuleType::EXCLUSIVE);
        bool good_input;

        if(ref_vm.count(m_requirement_pair.first))
        {
            good_input=(ref_vm.count(m_requirement_pair.first) && !ref_vm.count(m_requirement_pair.second));
        }

        else
        {
            good_input=true;
        }

        return good_input;
    }

    /*
    bool LaunchRule::_silent_parse(const po::variables_map &ref_vm) const
    {
        assert(m_type==RuleType::SILENT);
        bool good_input=false;

        return good_input;
    }

    bool LaunchRule::_required_parse(const po::variables_map &ref_vm) const
    {
        assert(m_type==RuleType::REQUIRED);
        bool good_input=false;

        return good_input;
    }
    */

    //****************************************************************************************************//
    
    /**
     * Creates two inclusion rules by permuting the order. If the given
     * target cannot go with the excluded option, the excluded option
     * cannot go with the target
     */

    void LaunchRuleList::add_exclusion(const std::string &target, const std::string &exclude)
    {
        m_strict_rule_table.push_back(LaunchRule(target, exclude, RuleType::EXCLUSIVE));
        m_strict_rule_table.push_back(LaunchRule(exclude, target, RuleType::EXCLUSIVE));

        return;
    }

    void LaunchRuleList::add_inclusion(const std::string &target, const std::string &include)
    {
        m_strict_rule_table.push_back(LaunchRule(target, include, RuleType::INCLUSIVE));

        return;
    }

    /**
     * Creates a list of inclusion rules, but only one of them must be obeyed. If none
     * of them are obeyed, or more than one are obeyed, then the ::parse routine will
     * return false.
     */

    void LaunchRuleList::add_any_inclusion(const std::string &target, const std::vector<std::string> &includes)
    {
        std::vector<LaunchRule> new_set;
        for(auto it=includes.begin(); it!=includes.end(); ++it)
        {
            new_set.push_back(LaunchRule(target, *it, RuleType::INCLUSIVE));
        }

        m_any_rule_table.push_back(new_set);

        return;
    }

    void LaunchRuleList::add_any_exclusion(const std::string &target, const std::vector<std::string> &excludes)
    {
        for(auto it=excludes.begin(); it!=excludes.end(); ++it)
        {
            add_exclusion(target, *it);
        }

        return;
    }
    

    /**
     * Runs through every possible combination of the given group of --suboptions
     * and creates an inclusion rule for each of them
     */

    //void LaunchRule::add_requirement(const std::vector<std::string> &codependent_group)
    //{
    //    for(auto it=codependent_group.begin(); it!=codependent_group.end(); ++it)
    //    {
    //        for(auto jt=codependent_group.begin(); jt!=codependent_group.end(); ++jt)
    //        {
    //        }
    //    }
    //    
    //    return;
    //}
    
    /**
     * Parse every rule and keep track of which ones were broken
     * by storing the indexes. Returns true only if every rule
     * passes.
     */

    bool LaunchRuleList::parse(const po::variables_map &ref_vm) const
    {
        bool all_rules_met=true;

        int index=0;
        for(auto it=m_strict_rule_table.begin(); it!=m_strict_rule_table.end(); ++it)
        {
            if(!it->parse(ref_vm))
            {
                all_rules_met=false;
                m_strict_broken_rules.push_back(index);
            }

            ++index;
        }

        index=0;
        for(auto it=m_any_rule_table.begin(); it!=m_any_rule_table.end(); ++it)
        {
            bool any_group_satisfied=false;
            for(auto jt=it->begin(); jt!=it->end(); ++jt)
            {
                if(jt->parse(ref_vm))
                {
                    any_group_satisfied=true;
                    break;
                }
            }

            if(!any_group_satisfied)
            {
                all_rules_met=false;
                m_any_broken_rules.push_back(index);
            }

            ++index;
        }


        return all_rules_met;
    }



}
