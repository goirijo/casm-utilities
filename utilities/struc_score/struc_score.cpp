#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace Utilities
{

void struc_score_initializer(po::options_description& struc_score_desc)
{
    UtilityProgramOptions::add_help_suboption(struc_score_desc);
    UtilityProgramOptions::add_output_suboption(struc_score_desc);

    struc_score_desc.add_options()("reference,r", po::value<fs::path>()->required(),
                    "POS.vasp like file to use as reference structure.");
    struc_score_desc.add_options()("mappable,m", po::value<std::vector<fs::path>>()->multitoken(),
                    "POS.vasp like file(s) to map, get structure score for."); 
    struc_score_desc.add_options()("batch,b", po::value<fs::path>(),
                    "Batch file containing list of structures files to get structure score for.");    
    struc_score_desc.add_options()("weight,w", po::value<double>()->default_value(0.5),
                    "Weight w in structure score: w*lattice_score + (1-w)*basis_score."); 
    struc_score_desc.add_options()("verbose,v","Print all three scores instead of only weighted."); 
    return;
}
}

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler struc_score_launch(argc, argv, struc_score_initializer);

    if (struc_score_launch.count("help"))
    {
        std::cout << struc_score_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        struc_score_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto reference_path = struc_score_launch.fetch<fs::path>("reference"); 
    auto map_reference_struc = CASM::Structure(reference_path);
    auto weight = struc_score_launch.fetch<double>("weight"); 
    auto verbose = struc_score_launch.count("verbose");

    std::vector<fs::path> mappable_paths;
    if (struc_score_launch.count("mappable"))
    {
        mappable_paths = struc_score_launch.fetch<std::vector<fs::path>>("mappable");
    }

    if (struc_score_launch.count("batch"))
    {
        std::ifstream batch_file(struc_score_launch.fetch<fs::path>("batch").c_str());
        std::string line;
        while(std::getline(batch_file, line))
        {
            mappable_paths.push_back(fs::path(line));
        }
    }

    if (mappable_paths.empty())
    {
        std::cout << "You must provide at least one structure to map." << std::endl; 
        return 2;
    } 
    
    std::vector<Rewrap::Structure> mappable_strucs;
    for (auto& path : mappable_paths)
    {
        mappable_strucs.push_back(Rewrap::Structure(path));
    }

    auto all_scores = Simplicity::structure_score(map_reference_struc, mappable_strucs);    
   
    std::string out_string = "";
    if (verbose)
    {
        out_string += "Structure\tLattice\tBasis\tWeighted\n";
    }
    for (int i=0; i < all_scores.size(); i++)
    {
        auto& sub_scores = all_scores[i];
        out_string += mappable_paths[i].string() + "\t";
        if (verbose)
        {
            out_string += std::to_string(sub_scores[0]) + "\t";
            out_string += std::to_string(sub_scores[1]) + "\t";
        }
        out_string += std::to_string(sub_scores[2]) + "\n";
     }
    
    if (struc_score_launch.count("output"))
    {
        auto out_path = struc_score_launch.fetch<fs::path>("output");
        std::ofstream file_out(out_path.string());
        file_out << out_string;
        file_out.close();
    }
    else
    {
        std::cout << out_string;
    }

    return 0;
}

