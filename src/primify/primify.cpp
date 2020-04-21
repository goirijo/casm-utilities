#include "CLI/Option.hpp"
#include "CLI/Validators.hpp"
#include <CLI/CLI.hpp>
#include <casmutils/definitions.hpp>
#include <casmutils/handlers.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>

namespace utilities
{
// TODO: Move somewhere common. You can make a custom callback so that
// fs::path is recognized as PATH instad of TEXT
CLI::Option* add_output_suboption(CLI::App* app, fs::path* output_path)
{
    auto* opt = app->add_option("-o,--output", *output_path, "Target output file");
    return opt;
}

} // namespace utilities

int main(int argc, char* argv[])
{
    CLI::App app;

    utilities::fs::path out_path;
    CLI::Option* out_path_opt = utilities::add_output_suboption(&app, &out_path);

    casmutils::fs::path super_path;
    CLI::Option* super_path_opt = app.add_option("-s,--superstructure",
                                                 super_path,
                                                 "POS.vasp like file that you want to get the primitive structure for.")
                                      ->required();
    super_path_opt->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);

    auto super_struc = casmutils::xtal::Structure::from_poscar(super_path);
    auto prim_struc = casmutils::xtal::make_primitive(super_struc);

    if (out_path_opt->count())
    {
        casmutils::xtal::write_poscar(prim_struc, out_path);
    }

    else
    {
        casmutils::xtal::print_poscar(prim_struc, std::cout);
    }

    return 0;
}
