#include <casmutils/definitions.hpp>
#include <casmutils/exceptions.hpp>
#include <casmutils/xtal/frankenstein.hpp>

#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/crystallography/io/VaspIO.hh>
#include <fstream>
#include <set>

namespace
{
using namespace casmutils;
/// This function alters the coordinates of the given struc to have fractional
/// coordinates between 0 and 1 only
void bring_coords_within(xtal::Structure* struc)
{
    throw except::NotImplemented();
    /* for (auto& site : struc->basis) */
    /* { */
    /*     site.within(); */
    /* } */
    return;
}
} // namespace

namespace casmutils
{
namespace frankenstein
{
void shift_coords_by(xtal::Structure* struc, const Eigen::Vector3d& shift_val)
{
    throw except::NotImplemented();
    /* for (auto& item : struc->basis) */
    /* { */
    /*     item += CASM::Coordinate(shift_val, struc->lattice(), CASM::FRAC); */
    /* } */
    /* bring_coords_within(struc); */
    /* return; */
}

std::pair<xtal::Structure, xtal::Structure> slice(const xtal::Structure& big_struc, double slice_loc, double tol)
{
    throw except::NotImplemented();
    /* // Copy the input, we'll do surgery on this guy */
    /* xtal::Structure cpy_big = big_struc; */
    /* bring_coords_within(&cpy_big); */

    /* // Create lattice for the bottom half of the structure */
    /* Eigen::Matrix3d lat_mat = big_struc.lattice().lat_column_mat(); */
    /* lat_mat.col(2) = lat_mat.col(2) * (slice_loc); */
    /* CASM::xtal::Lattice bottom_lat(lat_mat); */
    /* CASM::xtal::Structure bottom_struc(bottom_lat); */

    /* // Create lattice for the top half of the structure */
    /* lat_mat = big_struc.lattice().lat_column_mat(); */
    /* lat_mat.col(2) = lat_mat.col(2) * (1 - slice_loc); */
    /* CASM::xtal::Lattice top_lat(lat_mat); */
    /* CASM::xtal::Structure top_struc(top_lat); */

    /* // Loop over the basis of the input structure, and assign each */
    /* // site to either the top or bottom structures, shifting the */
    /* // sites accordingly */
    /* for (const auto& item : cpy_big.basis) */
    /* { */
    /*     // only move basis sites below slice pivot */
    /*     if (item.const_frac()(2) >= 0 && item.const_frac()(2) < slice_loc + tol) */
    /*     { */
    /*         auto coord = CASM::xtal::Coordinate(item.const_cart(), bottom_lat, CASM::CART); */
    /*         auto site = CASM::xtal::Site(coord, item.occ_name()); */
    /*         bottom_struc.basis.push_back(site); */
    /*     } */
    /*     else */
    /*     { */
    /*         CASM::Site new_site = item; */
    /*         // adjust c coord by slice location */
    /*         Eigen::Vector3d altered = new_site.const_frac(); */
    /*         altered(2) = new_site.const_frac()(2) - slice_loc; */
    /*         new_site.frac() = altered; */
    /*         auto coord = CASM::Coordinate(new_site.const_cart(), top_lat, CASM::CART); */
    /*         auto site = CASM::Site(coord, item.occ_name()); */
    /*         top_struc.basis.push_back(site); */
    /*     } */
    /* } */

    /* auto bottom_top = std::make_pair(xtal::Structure(bottom_struc), xtal::Structure(top_struc)); */
    /* bring_coords_within(&bottom_top.first); */
    /* bring_coords_within(&bottom_top.second); */
    /* return bottom_top; */
}

std::vector<xtal::Structure> _multi_slice(const xtal::Structure& big_struc, const Eigen::VectorXd& slice_locs,
                                          double tol)
{
    // Begin by performing the first slice
    auto struc_pair = slice(big_struc, slice_locs(0), tol);
    std::vector<xtal::Structure> slices;
    slices.push_back(struc_pair.first);

    // If that's the only slice left, you're done. Give the two
    // structures back
    if (slice_locs.size() == 1)
    {
        slices.push_back(struc_pair.second);
        return slices;
    }

    // You're not done, perform n-1 slices on the top structure
    Eigen::VectorXd subtract_from_locsations = Eigen::VectorXd::Constant(slice_locs.size() - 1, slice_locs(0));
    auto subtracted_locsations = slice_locs.tail(slice_locs.size() - 1) - subtract_from_locsations;

    // After subtracting the initial slice location, you also need to scale the slice locations
    // so the top structre slices are consistent with the initial structure slices
    auto after_piv = _multi_slice(struc_pair.second, subtracted_locsations / (1.0 - slice_locs(0)), tol);
    slices.insert(slices.end(), after_piv.begin(), after_piv.end());

    return slices;
}

std::vector<xtal::Structure> multi_slice(const xtal::Structure& big_struc, const std::set<double>& slice_locs,
                                         double tol)
{
    Eigen::VectorXd sanitized_slice_locs(slice_locs.size());

    int i = 0;
    for (const auto& slice : slice_locs)
    {
        if (slice < 0.0 || slice > 1.0)
        {
            throw except::UserInputMangle("Slice locations for frankenstein structures must be between 0.0 and 1.0");
        }

        sanitized_slice_locs(i) = slice;
        ++i;
    }

    return _multi_slice(big_struc, sanitized_slice_locs, tol);
}

std::vector<xtal::Structure> uniformly_slice(const xtal::Structure& big_struc, int n_pieces)
{
    std::set<double> slice_locs;
    for (int i = 1; i < n_pieces; i++)
    {
        slice_locs.insert((double)i / (double)n_pieces);
    }
    return multi_slice(big_struc, slice_locs, CASM::TOL);
}

xtal::Structure stack(const std::vector<xtal::Structure>& sub_strucs)
{
    throw except::NotImplemented();
    /* // Create a new lattice that has the same ab vectors. but summed up */
    /* // the c vectors of every structure */
    /* Eigen::Matrix3d stacked_lat_mat = sub_strucs[0].lattice().lat_column_mat(); */
    /* for (int i = 1; i < sub_strucs.size(); i++) */
    /* { */
    /*     Eigen::Matrix3d lat_mat = sub_strucs[i].lattice().lat_column_mat(); */
    /*     stacked_lat_mat.col(2) = stacked_lat_mat.col(2) + lat_mat.col(2); */
    /* } */
    /* CASM::Lattice stacked_lat(stacked_lat_mat); */

    /* // We now have a template structure with the right shape, we'll put the */
    /* // sites inside in a second. It already has the basis for the bottom of */
    /* // the stack */
    /* CASM::Structure stacked_struc(sub_strucs[0]); */
    /* stacked_struc.set_lattice(stacked_lat, CASM::CART); */

    /* // For each structure we stack, we'll take the basis, shift it up by the */
    /* // approprate amount, and stick it into our template stacked structure */
    /* Eigen::Vector3d c_shift = Eigen::Vector3d::Zero(); */
    /* for (int i = 1; i < sub_strucs.size(); i++) */
    /* { */
    /*     // determine appropriate c-axis shift for position in stacking */
    /*     c_shift = c_shift + sub_strucs[i].lattice().lat_column_mat().col(2); */
    /*     xtal::Structure cpy_i = sub_strucs[i]; */
    /*     bring_coords_within(&cpy_i); */

    /*     // Shift each site of the basis by the appropriate c shift, */
    /*     // and adds them to the stacked structure */
    /*     for (auto new_site : cpy_i.basis) */
    /*     { */
    /*         new_site.cart() += c_shift; */
    /*         new_site.set_lattice(stacked_lat, CASM::CART); */
    /*         stacked_struc.basis.push_back(new_site); */
    /*     } */
    /* } */

    /* auto rw_struc = xtal::Structure(stacked_struc); */
    /* bring_coords_within(&rw_struc); */
    /* return rw_struc; */
}

xtal::Structure vacuum_pack(const xtal::Structure& big_struc, std::array<bool, 3>& dirs, double padding)
{
    throw except::NotImplemented();
    /* xtal::Structure cpy_big = big_struc; */
    /* bring_coords_within(&cpy_big); */

    /* // We start by setting the boundaries at the maximum possible limits */
    /* // There is one pair for each lattice vector */
    /* std::vector<std::pair<double, double>> limits(3, std::make_pair(0.0, 1.0)); */

    /* // Going through each basis site, we determine the maximum and minimun boundaries, */
    /* // which correspond to the largest and smallest fractional coordinates of each lattice */
    /* // vector direction */
    /* for (auto& site : cpy_big.basis) */
    /* { */
    /*     auto coord = site.const_frac(); */
    /*     for (int i = 0; i < 3; i++) */
    /*     { */
    /*         if (coord(i) > limits[i].first) */
    /*         { */
    /*             limits[i].first = coord(i); */
    /*         } */
    /*         if (coord(i) < limits[i].second) */
    /*         { */
    /*             limits[i].second = coord(i); */
    /*         } */
    /*     } */
    /* } */

    /* // Shift the atoms so that the vacuum packed atoms are located at the origin */
    /* Eigen::Vector3d shift(-limits[0].second, -limits[1].second, -limits[2].second); */
    /* shift_coords_by(&cpy_big, shift); */

    /* // Reduce the lattice vectors so that the lattice only just encloses the atoms */
    /* Eigen::Matrix3d lat_mat = cpy_big.lattice().lat_column_mat(); */
    /* for (int i = 0; i < 3; i++) */
    /* { */
    /*     if (dirs[i]) */
    /*     { */
    /*         lat_mat.col(i) = lat_mat.col(i) * (limits[i].first - limits[i].second) + */
    /*                          padding * lat_mat.col(i) / lat_mat.col(i).norm(); */
    /*     } */
    /* } */

    /* // Set the lattice and you're done */
    /* xtal::Lattice lat(lat_mat); */
    /* cpy_big.set_lattice(lat, xtal::CART); */
    /* return cpy_big; */
}

xtal::Structure inflate(const xtal::Structure& struc, const std::array<double, 3>& padding)
{
    throw except::NotImplemented();
    // xtal::CasmStructure cpy_struc = struc;
    // Eigen::Matrix3d lat_mat = cpy_struc.lattice().lat_column_mat();

    //// Add padding to each lattice vector
    // for (int i = 0; i < 3; i++)
    //{
    //    lat_mat.col(i) = lat_mat.col(i) * (1.0 + padding[i] / lat_mat.col(i).norm());
    //}

    // cpy_struc.set_lattice(CASM::Lattice(lat_mat), CASM::CART);
    // xtal::Structure rw_struc(cpy_struc);
    // bring_coords_within(&rw_struc);
    // return rw_struc;
}

} // namespace frankenstein
} // namespace casmutils
