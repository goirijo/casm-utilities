
#include "lib/completer/handlers.hpp"
#include "lib/definitions.hpp"
#include "lib/exceptions.hpp"
#include "lib/launch/launchers.hpp"

#include <casm/CASM_global_definitions.hh>
#include <casm/app/AppIO.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/casm_io/jsonParser.hh>
#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/completer/Handlers.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/external/Eigen/Core>

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace casmUtilities
{
namespace bedazzleImpl
{

typedef std::pair<CASM::Structure, std::string> StrucStampinfoPair;
typedef std::map<std::string, StrucStampinfoPair> NameStrucMap;

void simple_structure_print(std::ostream &stream, const CASM::Structure &struc)
{
    CASM::VaspIO::PrintPOSCAR struc_printer(struc);
    struc_printer.sort();
    struc_printer.print(stream);
    return;
}

CASM::Structure stamped_structure(const CASM::SiteCluster &stamp, const CASM::Structure &background_struc)
{
    auto stamped_struc = background_struc;
    for (CASM::Index i = 0; i < stamp.size(); ++i)
    {
        for (CASM::Index j = 0; j < background_struc.basis.size(); ++j)
        {
            if (stamp[i] == background_struc.basis[j])
            {
                std::cout << stamp[i] << "        " << stamp[i].occ_name() << std::endl;
                std::cout << stamped_struc.basis[i] << "        " << stamped_struc.basis[i].occ_name() << std::endl;
                stamped_struc.basis[j] = stamp[i];
            }
        }
    }

    return stamped_struc;
}

std::vector<CASM::Index> linear_indexes_for_stamp(const CASM::Supercell &scel, const CASM::SiteCluster &stamp)
{
    std::vector<CASM::Index> linear_indexes;
    for (CASM::Index i = 0; i < stamp.size(); ++i)
    {
        linear_indexes.push_back(scel.get_linear_index(CASM::Coordinate(stamp[i]), CASM::TOL));
    }
    return linear_indexes;
}

std::string name_from_stamp(CASM::Index branch,
                            CASM::Index orbit,
                            CASM::Index decor,
                            const CASM::Array<CASM::Array<int>> &decor_map,
                            const CASM::SiteCluster &stamp)
{
    std::string name = "dec";
    name += std::to_string(branch) + ".";
    name += std::to_string(orbit) + ".";
    name += std::to_string(decor);

    if (decor_map[decor].size() > 0)
    {
        name += "_";
    }

    for (CASM::Index j = 0; j < decor_map[decor].size(); ++j)
    {
        name += stamp[j].allowed_occupants()[decor_map[decor][j]];
    }

    return name;
}

bool import_self_config(const CASM::Configuration &background_config, std::string &import_name)
{
    CASM::Index new_ix;
    CASM::Supercell::permute_const_iterator perm_it;
    bool new_config = background_config.get_supercell().add_config(background_config, new_ix, perm_it);
    import_name = background_config.get_supercell().get_config(new_ix).name();

    return new_config;
}

void print_stamp_info(std::ostream &stream,
                      const CASM::SiteCluster &stamp,
                      CASM::Index decor,
                      const CASM::Array<CASM::Array<int>> &decor_map)
{
    stream << "Points: " << stamp.size() << std::endl;
    stream << "Max length: " << stamp.max_length() << std::endl;
    stream << "Min length: " << stamp.min_length() << std::endl;

    for (int i = 0; i < stamp.size(); ++i)
    {
        stamp[i].print(stream);
        stream << " :: " << stamp[i].allowed_occupants()[decor_map[decor][i]] << std::endl;
    }
}

CASM::Configuration &stamp_configuration(CASM::Index decor,
                                         const std::vector<CASM::Index> &linear_indexes,
                                         const CASM::Array<CASM::Array<int>> &decor_map,
                                         CASM::Configuration &background_config)
{
    for (CASM::Index j = 0; j < decor_map[decor].size(); ++j)
    {
        background_config.set_occ(linear_indexes[j], decor_map[decor][j]);
    }
}

NameStrucMap &insert_structure(NameStrucMap &dec_to_configname,
                               CASM::Index branch,
                               CASM::Index orbit,
                               CASM::Index decor,
                               const CASM::Array<CASM::Array<int>> &decor_map,
                               const std::vector<CASM::Index> &linear_indexes,
                               CASM::Configuration background_config,
                               const CASM::SiteCluster &stamp)
{
    stamp_configuration(decor, linear_indexes, decor_map, background_config);

    std::string new_config_name;
    bool new_config = import_self_config(background_config, new_config_name);

    if (new_config)
    {
        auto new_dec_name = name_from_stamp(branch, orbit, decor, decor_map, stamp);
        std::stringstream stamp_desc_stream;
        print_stamp_info(stamp_desc_stream, stamp, decor, decor_map);

        StrucStampinfoPair info_pair(background_config.get_supercell().superstructure(background_config),
                                     stamp_desc_stream.str());
        dec_to_configname[new_dec_name] = info_pair;
    }

    /* background_config.set_occupation(background_occ); */

    return dec_to_configname;
}

NameStrucMap self_stamp_primclex(const CASM::SiteOrbitree &background_tree, CASM::Configuration &background_config)
{
    // Map the bedazzle names to the config names in the PrimClex
    NameStrucMap dec_to_struc;

    // Save the original occupation of the background to reset later
    const auto background_occ = background_config.occupation();

    // Naturally we are using the Array as a container, because we're better
    // than std::vector dammit!
    for (CASM::Index b = 0; b < background_tree.size(); ++b)
    {
        for (CASM::Index o = 0; o < background_tree[b].size(); ++o)
        {
            auto stamp = background_tree[b][o].prototype;
            auto linear_indexes = linear_indexes_for_stamp(background_config.get_supercell(), stamp);
            auto decor_map = stamp.get_full_decor_map();

            for (CASM::Index i = 0; i < decor_map.size(); ++i)
            {
                insert_structure(dec_to_struc, b, o, i, decor_map, linear_indexes, background_config, stamp);
            }
        }
    }

    return dec_to_struc;
}

void write_dec_dirs(const NameStrucMap &dec_to_struc,
                    CASM::fs::path write_dir,
                    std::string pos_filename,
                    std::string stamp_filename)
{
    for (auto &dec : dec_to_struc)
    {
        auto target_dir = write_dir / dec.first;
        auto target_poscar = target_dir / pos_filename;
        auto target_stamp = target_dir / stamp_filename;

        if (!CASM::fs::exists(target_dir))
        {
            CASM::fs::create_directories(target_dir);
        }

        else
        {
            throw OverwriteException(target_dir.string());
        }

        std::ofstream pos_stream, stamp_stream;
        pos_stream.open(target_poscar.string());
        simple_structure_print(pos_stream, dec.second.first);
        pos_stream.close();

        stamp_stream.open(target_stamp.string());
        stamp_stream << dec.second.second << std::endl;
        stamp_stream.close();
    }
    return;
}
}

std::string bedazzle_launcher_name() { return "bedazzle"; }

LaunchRuleList bedazzle_initializer(po::options_description &bedazzle_desc)
{
    utilityProgramOptions::add_help_suboption(bedazzle_desc);

    // clang-format off
    bedazzle_desc.add_options()
        ("cspecs,c", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path())->required(),
                                "bspecs.json like file for defining cluster defects.")
        ("prim,p", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path())->required(),
        "prim.json like file that defines the primitive cell that the background structure is relative to.")
        ("background-struc,b", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path())->required(),
        "POS.vasp like file that is the configuration that you want to decorate. Must be a superstructure of the "
        "provided primitive strucutre.")
        ("stamp-filename,s", po::value<std::string>()->default_value("stamp.txt"),
                                         "Filename for what will contain the information of the applied defects.")
        ("poscar-filename,f", po::value<std::string>()->default_value("POSCAR"),
        "Filename for the bedazzled structure files.")
        ("lattice-struc,l", po::value<fs::path>()->value_name(CASM::Completer::ArgHandler::path()),"POS.vasp like file that defines what the final lattice of the structures should be. If not given the lattice of --prim is kept.")
        ("target-dir,t", po::value<fs::path>()->default_value("./"),"Target directory to where the output should go");

    LaunchRuleList bedazzle_rules;

    return bedazzle_rules;
}

void bedazzle_utility_launch(int argc, char *argv[])
{
    Launcher bedazzle_launch(argc, argv, bedazzle_launcher_name(), bedazzle_initializer);

    if (bedazzle_launch.count("help"))
    {
        std::cout << bedazzle_launch.utility().desc() << std::endl;
        return;
    }

    try
    {
        bedazzle_launch.notify();
    }

    catch (boost::program_options::required_option &e)
    {
        std::cerr << e.what() << std::endl;
        return;
    }

    auto cspecs_file = bedazzle_launch.fetch<fs::path>("cspecs");
    auto prim_file = bedazzle_launch.fetch<fs::path>("prim");
    auto background_file = bedazzle_launch.fetch<fs::path>("background-struc");

    // First read the background structure, which might have been relaxed
    // and oriended some weird way, and a prim structure that specifies
    // how the occupants can behave

    CASM::BasicStructure<CASM::Site> background_struc;
    CASM::fs::ifstream background_file_stream(background_file);
    background_struc.read(background_file_stream);

    CASM::BasicStructure<CASM::Site> prim_struc(CASM::read_prim(prim_file));

    // Construct a primclex out of the primitive structure, this will
    // serve to standardize the structure we want to decorate
    CASM::Structure init_struc(prim_struc);
    CASM::PrimClex pclex(init_struc);

    // Default values for ConfigMapper
    double lat_weight = 0.5;
    CASM::ConfigMapper cmapper(pclex, lat_weight);

    std::string background_name;
    CASM::jsonParser tmp_json;
    std::vector<CASM::Index> tmp_Ix_vector;
    Eigen::Matrix3d tmp_3d;

    // Test importing structures
    bool mapped;
    mapped = cmapper.import_structure_occupation(background_struc, background_name, tmp_json, tmp_Ix_vector, tmp_3d);

    // Create an Oribtree
    CASM::fs::path bspecs_path("./perturbspecs.json");
    CASM::jsonParser bspecs(bspecs_path);
    // This is why I hate CASM. Use the structure to make a config to get the
    // supercell to make the structure from the config
    auto background_config = pclex.configuration(background_name);
    CASM::Structure background_struc_recast(background_config.get_supercell().superstructure(background_config));
    CASM::Structure prim_struc_recast(prim_struc);
    CASM::SiteOrbitree background_tree = make_orbitree(background_struc_recast, bspecs);

    background_tree.print(std::cout);

    // Let the bedazzle begin
    bedazzleImpl::NameStrucMap dec_to_struc = bedazzleImpl::self_stamp_primclex(background_tree, background_config);
    auto target_dir = bedazzle_launch.fetch<fs::path>("target-dir");
    auto stamp_filename = bedazzle_launch.fetch<std::string>("stamp-filename");
    auto poscar_filename = bedazzle_launch.fetch<std::string>("poscar-filename");

    //reset the lattice for all structures to the original
    if(bedazzle_launch.count("lattice-struc"))
    {
        auto lat_file = bedazzle_launch.fetch<fs::path>("lattice-struc");
        CASM::Structure lat_set_struc(lat_file);
        for (auto &dec : dec_to_struc)
        {
            dec.second.first.set_lattice(lat_set_struc.lattice(), CASM::FRAC);
        }

    }

    bedazzleImpl::write_dec_dirs(dec_to_struc, target_dir.string(),poscar_filename,stamp_filename);
    return;
}
}
