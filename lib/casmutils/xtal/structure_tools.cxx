#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/crystallography/BasicStructureTools.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/SuperlatticeEnumerator.hh>
#include <casm/crystallography/io/VaspIO.hh>
#include <casm/strain/StrainConverter.hh>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <fstream>

namespace
{
// surface area of a given lattice
double lattice_surface_area(const casmutils::xtal::Lattice& lat)
{
    Eigen::Vector3d a = lat[0];
    Eigen::Vector3d b = lat[1];
    Eigen::Vector3d c = lat[2];

    double ab = a.cross(b).norm();
    double bc = b.cross(c).norm();
    double ca = c.cross(b).norm();

    return std::abs(ab) + std::abs(bc) + std::abs(ca);
}

// score for determining level of boxiness, borrowed from John Goiri
double boxy_score(const casmutils::xtal::Lattice& lat)
{
    // Less surface area per volume means more boxy
    // i.e. more volume per surface area means more boxy
    return std::abs(lat.volume()) / lattice_surface_area(lat);
}
} // namespace

namespace casmutils
{
namespace xtal
{
Structure make_primitive(const Structure& input)
{
    const auto& casted_input = input.__get<CASM::xtal::BasicStructure>();
    Structure true_prim(CASM::xtal::make_primitive(casted_input));
    return true_prim;
}

Structure make_niggli(const Structure& non_niggli)
{
    throw except::NotImplemented();
    /* CasmStructure niggli = non_niggli; */
    /* CASM::Lattice lat_niggli = CASM::xtal::niggli(non_niggli.lattice(), CASM::TOL); */
    /* niggli.set_lattice(lat_niggli, CASM::CART); */
    /* return niggli; */
}

void make_niggli(Structure* non_niggli)
{
    throw except::NotImplemented();
    /* CASM::Lattice lat_niggli = CASM::xtal::niggli(non_niggli->lattice(), CASM::TOL); */
    /* non_niggli->set_lattice(lat_niggli, CASM::CART); */
    /* return; */
}

void print_poscar(const Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(printable.__get<CASM::xtal::SimpleStructure>());
    p.sort();
    p.print(outstream);
    return;
}

void write_poscar(const Structure& printable, const fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}

Structure make_super_structure(const Structure& struc, const Eigen::Matrix3i& col_transf_mat)
{
    auto lattice_mat = struc.lattice().column_vector_matrix();
    // had to cast the transformation matrix to double as Eigen does not allow mixing matrix types
    CASM::xtal::Lattice suplat(lattice_mat * col_transf_mat.cast<double>());
    return Structure(struc.__get<CASM::xtal::BasicStructure>().create_superstruc(suplat));
}

void apply_deformation(Structure* struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    Lattice strained_lattice(deformation_tensor * struc_ptr->lattice().column_vector_matrix());
    struc_ptr->set_lattice(strained_lattice, FRAC);
    return;
}

Structure apply_deformation(const Structure& struc_ptr, const Eigen::Matrix3d& deformation_tensor)
{
    Structure copy_struc(struc_ptr);
    apply_deformation(&copy_struc, deformation_tensor);
    return copy_struc;
}

void apply_strain(Structure* struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
{
    std::set<std::string> allowed_strain_metrics = {"GL", "B", "H", "EA"};
    if (allowed_strain_metrics.count(mode))
    {
        CASM::StrainConverter converter(mode);
        auto strain_tensor = converter.rollup_E(unrolled_strain);
        auto deformation_tensor = converter.strain_metric_to_F(strain_tensor);
        apply_deformation(struc_ptr, deformation_tensor);
    }
    else
    {
        throw except::UserInputMangle("Unrecognized mode. Allowed strain metrics modes are GL, B, H and EA");
    }
    return;
}

Structure apply_strain(const Structure& struc_ptr, const Eigen::VectorXd& unrolled_strain, const std::string& mode)
{
    Structure copy_struc(struc_ptr);
    apply_strain(&copy_struc, unrolled_strain, mode);
    return copy_struc;
}

std::vector<std::pair<double, double>> structure_score(const Structure& map_reference_struc,
                                                       const std::vector<Structure>& mappable_struc_vec)
{
    throw except::NotImplemented();
    /* for (const auto& struc : mappable_struc_vec) */
    /* { */
    /*     if (struc.basis.size() != map_reference_struc.basis.size()) */
    /*     { */
    /*         throw except::BasisMismatch(); */
    /*     } */
    /* } */

    /* // get prim and make PrimClex */
    /* auto ref_prim = CASM::Structure(make_primitive(map_reference_struc)); */
    /* auto pclex = extend::quiet_primclex(ref_prim); */

    /* // mapping setup */
    /* int options = 2;      // robust mapping */
    /* double vol_tol = 0.5; // not used */
    /* double weight = 0.5;  // not used */
    /* CASM::ConfigMapper configmapper(pclex, weight, vol_tol, options, CASM::TOL); */

    /* CASM::jsonParser out;          // mapping output */
    /* std::string name;              // not used */
    /* std::vector<CASM::Index> best; // not used */
    /* Eigen::Matrix3d cart_op;       // not used */
    /* bool update = false; */

    /* std::vector<std::pair<double, double>> all_scores; */
    /* for (const auto& struc : mappable_struc_vec) */
    /* { */
    /*     // map it */
    /*     Structure mappable_copy(struc); // can't be const, make copy */
    /*     configmapper.import_structure_occupation(mappable_copy, name, out, best, cart_op, update); */
    /*     double basis = out["best_mapping"]["basis_deformation"].get<double>(); */
    /*     double lattice = out["best_mapping"]["lattice_deformation"].get<double>(); */
    /*     all_scores.emplace_back(lattice, basis); */
    /* } */

    /* return all_scores; */
}

std::pair<double, double> structure_score(const Structure& map_reference_struc, const Structure& mappable_struc)
{
    std::vector<Structure> one = {mappable_struc};
    // just calls vector version of function
    return structure_score(map_reference_struc, one).back();
}

// Finds the superstructure with the highest volume/surface_area
// Assuming that the input has structures of same volume
std::vector<Structure>::size_type boxiest_structure_index(const std::vector<Structure>& candidate_structures)
{
    // TODO: throw exception on empty vector
    double running_score = 0;
    std::vector<Structure>::size_type ix = 0;
    std::vector<Structure>::size_type best_ix = ix;
    for (const auto& scel : candidate_structures)
    {
        double candidate_score = boxy_score(scel.lattice());
        if (candidate_score > running_score)
        {
            running_score = candidate_score;
            best_ix = ix;
        }
        ++ix;
    }
    return best_ix;
}

// Find the boxiest superstructure per volume for range of volumes
Structure make_boxiest_superstructure_of_volume(const Structure& structure, const int volume)
{
    std::vector<Structure> same_vol_scels = make_superstructures_of_volume(structure, volume);
    return same_vol_scels[boxiest_structure_index(same_vol_scels)];
}

std::vector<Structure> make_superstructures_of_volume(const Structure& structure, const int volume)
{
    std::vector<Structure> all_superstructures;
    /* CASM::xtal::ScelEnumProps enum_props(volume, volume+1); */
    /* CASM::xtal::SuperlatticeEnumerator lat_enumerator(structure.lattice(), enum_props, CASM::TOL); */

    /* for (const auto& lat : lat_enumerator) */
    /* { */
    /*     Structure super = structure.create_superstruc(lat); */
    /*     simplicity::make_niggli(&super); */
    /*     all_superstructures.emplace_back(std::move(super)); */
    /* } */

    throw except::NotImplemented();

    return all_superstructures;
}
} // namespace xtal
} // namespace casmutils
