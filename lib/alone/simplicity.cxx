#include <iostream>

#include <casm/app/AppIO.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/Structure.hh>

namespace simple
{
void structure_print(std::ostream &stream, const CASM::Structure &struc)
{
    CASM::VaspIO::PrintPOSCAR struc_printer(struc);
    struc_printer.sort();
    struc_printer.print(stream);
    return;
}

CASM::PrimClex primclex_from_path(CASM::fs::path prim_path)
{
    CASM::BasicStructure<CASM::Site> prim_struc(CASM::read_prim(prim_path));
    CASM::Structure init_struc(prim_struc);
    CASM::PrimClex pclex(init_struc);

    return pclex;
}

CASM::BasicStructure<CASM::Site> basic_structure_from_path(CASM::fs::path struc_path)
{
    CASM::fs::ifstream base_struc_stream(struc_path);
    CASM::BasicStructure<CASM::Site> base_struc;
    base_struc.read(base_struc_stream);

    return base_struc;
}

bool primclex_import_basic_structure(CASM::PrimClex &pclex, CASM::BasicStructure<CASM::Site> new_struc)
{
    double lat_weight = 0.5;
    CASM::ConfigMapper cmapper(pclex, lat_weight);
    std::string struc_name;
    CASM::jsonParser tmp_json;
    std::vector<CASM::Index> tmp_Ix_vector;
    Eigen::Matrix3d tmp_3d;

    bool mapped = cmapper.import_structure_occupation(new_struc, struc_name, tmp_json, tmp_Ix_vector, tmp_3d);
    return mapped;
}

Eigen::MatrixXi easy_matrix_round(const Eigen::MatrixXd &mat)
{
    Eigen::MatrixXi int_mat(mat.rows(), mat.cols());
    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            int_mat(i, j) = int(mat(i, j) + 0.5);
        }
    }
    return int_mat;
}
}
