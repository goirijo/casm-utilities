#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

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
    fs::path::string_type::size_type max_path_length=0;
    for (auto& path : mappable_paths)
    {
        mappable_strucs.push_back(Rewrap::Structure(path));
        max_path_length=std::max(max_path_length,path.size());
    }

    auto all_scores = Simplicity::structure_score(map_reference_struc, mappable_strucs);    
    
    std::ostream* out_stream_ptr=&std::cout;
    std::ofstream specified_out_stream;
    if (struc_score_launch.count("output"))
    {
        auto out_path = struc_score_launch.fetch<fs::path>("output");
        specified_out_stream.open(out_path.c_str());
        out_stream_ptr=&specified_out_stream;
    }
    auto& out_stream=*out_stream_ptr;

    out_stream<<std::left<<std::setw(max_path_length+4)<<"Structure";
    out_stream<<std::left<<std::setw(16)<<"Lattice";
    out_stream<<std::left<<std::setw(16)<<"Basis";
    out_stream<<std::left<<std::setw(16)<<"Weighted"<<std::endl;

    assert(all_scores.size()==mappable_paths.size());
    int path_ix=0;
    for (const auto& lat_basis_score : all_scores)
    {
        auto lat_score=lat_basis_score.first;
        auto basis_score=lat_basis_score.second;
        auto weighted_score=weight*lat_score+(1-weight)*basis_score;

        out_stream<<std::left<<std::setw(max_path_length+4)<<mappable_paths[path_ix].string();
        out_stream<<std::left<<std::setw(16)<<std::setprecision(8)<<lat_score;
        out_stream<<std::left<<std::setw(16)<<std::setprecision(8)<<basis_score;
        out_stream<<std::left<<std::setw(16)<<std::setprecision(8)<<weighted_score<<std::endl;

        ++path_ix;
     }

    //why though
    if (struc_score_launch.count("output"))
    {
        specified_out_stream.close();
    }

    return 0;
}

