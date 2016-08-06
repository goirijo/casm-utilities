#include <iterator>

#include "lib/launch/rules.hpp"

namespace casmUtilities
{
    /**
     * Creates two inclusion rules by permuting the order. If the given
     * target cannot go with the excluded option, the excluded option
     * cannot go with the target
     */

    void LaunchRule::add_exclusion(std::string &target, std::string &exclude)
    {
        m_exclusion_table.push_back(std::make_pair(target, exclude));
        m_exclusion_table.push_back(std::make_pair(exclude, target));

        return;
    }

    void LaunchRule::add_inclusion(std::string &target, std::string &include)
    {
        m_inclusion_table.push_back(std::make_pair(target, include));

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
}
