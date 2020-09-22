#include "casmutils/xtal/coordinate.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <casmutils/mush/slab.hpp>
#include <casmutils/mush/twist.hpp>
#include <casmutils/xtal/structure.hpp>
/* #include "multishift/slab.hpp" */
#include <casm/crystallography/SuperlatticeEnumerator.hh>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <cassert>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
using namespace casmutils;

Eigen::Matrix2d make_2d_column_matrix(const xtal::Lattice lat3d)
{
    assert(CASM::almost_equal(lat3d.column_vector_matrix()(2, 0), 0.0, 1e-10));
    assert(CASM::almost_equal(lat3d.column_vector_matrix()(2, 1), 0.0, 1e-10));

    return lat3d.column_vector_matrix().block<2, 2>(0, 0);

    Eigen::Vector2d a2d(lat3d.a()(0), lat3d.a()(1));
    Eigen::Vector2d b2d(lat3d.b()(0), lat3d.b()(1));

    Eigen::Matrix2d lat_mat2d;
    lat_mat2d.col(0) = a2d;
    lat_mat2d.col(1) = b2d;

    return lat_mat2d;
};

Eigen::Matrix3d make_3d_column_matrix(const Eigen::Matrix2d& lat2d)
{
    Eigen::Matrix3d lat3d = Eigen::Matrix3d::Identity();
    lat3d.block<2, 2>(0, 0) = lat2d;
    return lat3d;
}

double coarse_round(double x)
{
    double trucated = std::round(100000.0 * x) / 100000.0;
    return std::round(trucated);
}

bool is_within_voronoi(const Eigen::Vector3d& v, const xtal::Lattice& lat)
{
    xtal::Coordinate vw(v);
    vw.bring_within_wigner_seitz(lat);
    return almost_equal(v, vw.cart(), 1e-13);
}

} // namespace

namespace casmutils
{
namespace mush
{
Eigen::Matrix3d make_twist_rotation_matrix(const xtal::Lattice& lat, double degrees)
{
    double rad = M_PI * degrees / 180.0;
    Eigen::Matrix3d z_rot_mat;
    z_rot_mat << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;

    Eigen::Matrix3d standard_to_lat_span = slab_unit_vectors(lat);
    Eigen::Matrix3d lat_span_to_standard = standard_to_lat_span.inverse();

    return standard_to_lat_span * z_rot_mat * lat_span_to_standard;
}

xtal::Lattice make_twisted_lattice(const xtal::Lattice& lat, double degrees)
{
    Eigen::Matrix3d twist_matrix = make_twist_rotation_matrix(lat, degrees);
    return make_transformed_lattice(lat, twist_matrix);
}

std::pair<Eigen::Matrix3l, Eigen::Matrix3d> approximate_integer_transformation(const xtal::Lattice& L,
                                                                               const xtal::Lattice& M)
{
    Eigen::Matrix3d Td = L.column_vector_matrix().inverse() * M.column_vector_matrix();

    Eigen::Matrix3l T = CASM::lround(Td);
    Eigen::Matrix3d E = Td - T.cast<double>();

    return std::make_pair(T, E);
}

MoireLattice::MoireLattice(const xtal::Lattice& lat, double degrees)
    : input_lattice(lat),
      input_degrees(degrees),
      aligned_lattice(make_prismatic_lattice(make_aligned_lattice(input_lattice))),
      reciprocal_aligned_lattice(xtal::make_reciprocal(aligned_lattice)),
      rotated_lattice(make_twisted_lattice(aligned_lattice, input_degrees)),
      reciprocal_rotated_lattice(xtal::make_reciprocal(rotated_lattice)),
      full_reciprocal_difference(this->calculate_reciprocal_difference()),
      aligned_brillouin_zone_reciprocal_difference(
          this->bring_vectors_into_voronoi(full_reciprocal_difference, reciprocal_aligned_lattice)),
      rotated_brillouin_zone_reciprocal_difference(
          this->bring_vectors_into_voronoi(full_reciprocal_difference, reciprocal_rotated_lattice)),
      aligned_moire_lattice(make_moire_lattice_from_reciprocal_difference(aligned_brillouin_zone_reciprocal_difference,
                                                                          aligned_lattice.c())),
      rotated_moire_lattice(make_moire_lattice_from_reciprocal_difference(rotated_brillouin_zone_reciprocal_difference,
                                                                          aligned_lattice.c()))
{
    // calculate overlap
    for (int i = 0; i < 2; ++i)
    {
        this->is_within_brillouin_zone_overlap[LATTICE::ALIGNED][i] =
            this->is_within_voronoi(aligned_brillouin_zone_reciprocal_difference.col(i), reciprocal_rotated_lattice);
        this->is_within_brillouin_zone_overlap[LATTICE::ROTATED][i] =
            this->is_within_voronoi(rotated_brillouin_zone_reciprocal_difference.col(i), reciprocal_aligned_lattice);
    }
}

bool MoireLattice::is_within_voronoi(const Eigen::Vector2d& v, const xtal::Lattice& lat) const
{
    Eigen::Vector3d v3(v(0), v(1), 0);
    return ::is_within_voronoi(v3, lat);
}

Eigen::Matrix2d MoireLattice::calculate_reciprocal_difference() const
{
    Eigen::Matrix3d diff = this->reciprocal_rotated_lattice.column_vector_matrix() -
                           this->reciprocal_aligned_lattice.column_vector_matrix();

    Eigen::Vector3d z3 = Eigen::Vector3d::Zero();

    // There should be no difference in the c vector
    assert(almost_equal(diff.col(2), z3));
    // There should be no z components for anything
    assert(almost_equal(diff.row(2), z3.transpose()));

    return diff.block<2, 2>(0, 0);
}

Eigen::Matrix2d MoireLattice::bring_vectors_into_voronoi(const Eigen::Matrix2d& col_vectors, const xtal::Lattice& lat)
{
    // You must provide a prismatic aligned lattice for this to word
    if (!CASM::almost_equal(lat.a()(2), 0.0) && !CASM::almost_equal(lat.b()(2), 0.0) &&
        !CASM::almost_equal(lat.c()(0), 0.0) && !CASM::almost_equal(lat.c()(1), 0.0))
    {
        throw std::runtime_error(
            "The provided lattice must be prismatic (perpendicular c vector) and aligned in the xy plane");
    }

    // TODO: Bring within voronoi. Function is missing in cu
    xtal::Coordinate ka(col_vectors(0, 0), col_vectors(1, 0), 0);
    xtal::Coordinate kb(col_vectors(0, 1), col_vectors(1, 1), 0);

    ka.bring_within_wigner_seitz(lat);
    kb.bring_within_wigner_seitz(lat);

    Eigen::Matrix2d col_vectors_within;
    col_vectors_within.col(0) = ka.cart().head(2);
    col_vectors_within.col(1) = kb.cart().head(2);

    return col_vectors_within;
}

xtal::Lattice MoireLattice::make_right_handed_by_ab_swap(const xtal::Lattice& lat)
{
    if (lat.column_vector_matrix().determinant() < 0)
    {
        return xtal::Lattice(lat.b(), lat.a(), lat.c());
    }
    return lat;
}

xtal::Lattice MoireLattice::make_moire_lattice_from_reciprocal_difference(const Eigen::Matrix2d diff,
                                                                          const Eigen::Vector3d& real_c_vector)
{
    Eigen::Matrix3d phony_recip_mat = Eigen::Matrix3d::Identity();
    phony_recip_mat.block<2, 2>(0, 0) = diff;
    xtal::Lattice phony_recip_lat(phony_recip_mat.col(0), phony_recip_mat.col(1), phony_recip_mat.col(2));

    xtal::Lattice phony_real_moire_lattice = xtal::make_reciprocal(phony_recip_lat);
    Eigen::Matrix3d final_real_moire_mat = phony_real_moire_lattice.column_vector_matrix();
    final_real_moire_mat.col(2) = real_c_vector;

    xtal::Lattice moire(final_real_moire_mat.col(0), final_real_moire_mat.col(1), final_real_moire_mat.col(2));
    return make_right_handed_by_ab_swap(moire);
}

MoireApproximant::MoireApproximant(const xtal::Lattice& moire_lat,
                                   const xtal::Lattice& aligned_lat,
                                   const xtal::Lattice& rotated_lat)
    : approximate_moire_lattice(this->default_lattice())
{
    std::unordered_map<LATTICE, const xtal::Lattice*> real_lattices;
    real_lattices.emplace(LATTICE::ALIGNED, &aligned_lat);
    real_lattices.emplace(LATTICE::ROTATED, &rotated_lat);

    // Figure out the integer transformations
    for (LATTICE lat : {LATTICE::ALIGNED, LATTICE::ROTATED})
    {
        const auto& M = moire_lat;
        const auto& L = *real_lattices[lat];

        std::tie(approximate_moire_integer_transformations[lat], approximate_moire_integer_transformation_errors[lat]) =
            approximate_integer_transformation(L, M);
    }

    // Find the approximate Moire lattice
    auto aligned_S =
        xtal::make_superlattice(aligned_lat, approximate_moire_integer_transformations[LATTICE::ALIGNED].cast<int>());
    auto rotated_S =
        xtal::make_superlattice(rotated_lat, approximate_moire_integer_transformations[LATTICE::ROTATED].cast<int>());
    Eigen::Matrix3d S_bar = (aligned_S.column_vector_matrix() + rotated_S.column_vector_matrix()) / 2.0;
    this->approximate_moire_lattice = xtal::Lattice(S_bar);

    // Determine the strain involved to make things purrfect
    Eigen::Matrix3d aligned_F = S_bar * aligned_S.column_vector_matrix().inverse();
    Eigen::Matrix3d rotated_F = S_bar * rotated_S.column_vector_matrix().inverse();
    approximation_deformations[LATTICE::ALIGNED] = aligned_F;
    approximation_deformations[LATTICE::ROTATED] = rotated_F;

    // Determine the deformed tiling units
    for (LATTICE lat : {LATTICE::ALIGNED, LATTICE::ROTATED})
    {
        const Eigen::Matrix3d& L = real_lattices[lat]->column_vector_matrix();
        Eigen::Matrix3d L_prime = S_bar * approximate_moire_integer_transformations[lat].cast<double>().inverse();
        approximate_lattices.emplace(lat, L_prime);
    }
}

xtal::Lattice make_prismatic_lattice(const xtal::Lattice& lat)
{
    Eigen::Vector3d orthogonal_unit = lat.a().cross(lat.b()).normalized();
    Eigen::Vector3d new_c = lat.c().dot(orthogonal_unit) * orthogonal_unit;
    return xtal::Lattice(lat.a(), lat.b(), new_c);
}

DeformationReport::DeformationReport(const Eigen::Matrix3d& _deformation):deformation(_deformation)
{
    std::tie(rotation,strain)=xtal::polar_decomposition(deformation);
    this->rotation_angle=std::atan2(rotation(1,0),rotation(0,0))*180.0/M_PI;

    Eigen::Matrix3d E=strain-Eigen::Matrix3d::Identity();
    this->strain_metrics[0]=(1.0/std::sqrt(2))*(E(0,0)+E(1,1));
    this->strain_metrics[1]=(1.0/std::sqrt(2))*(E(0,0)-E(1,1));
    this->strain_metrics[2]=std::sqrt(2)*E(0,1);

    const auto& eta=strain_metrics;
    dilation_strain=eta[0];
    deviatoric_strain=std::sqrt(eta[1]*eta[1]+eta[2]*eta[2]);

    if(!almost_zero(E.col(2)) || !almost_zero(E.row(2)))
    {
        throw std::runtime_error("Deformation matrix extends beyond the xy subspace");
    }

}

MoireLatticeReport MoireGenerator::generate(ZONE bz, LATTICE layer) const
{
    const Eigen::Matrix3l& true_moire_scel_mat = (bz == ZONE::ALIGNED) ? transformation_matrix_to_super_aligned_moire
                                                            : transformation_matrix_to_super_rotated_moire;

    const auto& approximant = (bz == ZONE::ALIGNED) ? aligned_moire_approximant : rotated_moire_approximant;
    /* return MoireLatticeReport(bz, */
    /*                           layer, */
    /*                           /1* angle, *1/ */
    /*                           moire.moire(bz), */
    /*                           true_moire_scel_mat, */
    /*                           approximant.approximate_moire_lattice, */
    /*                           approximant.approximate_lattices.at(layer), */
    /*                           approximant.approximate_moire_integer_transformations.at(layer), */
    /*                           Eigen::Matrix3l(), */
    /*                           Eigen::Matrix3d()); */
    return MoireLatticeReport(bz,
                              layer,
                              /* angle, */
                              moire.moire(bz),
                              true_moire_scel_mat,
                              approximant.approximate_moire_lattice,
                              approximant.approximate_lattices.at(layer),
                              approximant.approximate_moire_integer_transformations.at(layer),
                              approximant.approximate_moire_integer_transformation_errors.at(layer),
                              approximant.approximation_deformations.at(layer));
}

xtal::Lattice MoireGenerator::make_reduced_cell(const xtal::Lattice& lat)
{
    auto raw_reduced = lat.__get().reduced_cell2();

    Eigen::Matrix3d T;
    // turns a,b,c to b,c,a
    Eigen::Matrix3d P;
    P << 0, 1, 0, 0, 0, 1, 1, 0, 0;
    // Swap a and b
    Eigen::Matrix3d S;
    S << 0, 1, 0, 1, 0, 0, 0, 0, 1;

    if (almost_equal(lat.c(), raw_reduced[2]))
    {
        T = Eigen::Matrix3d::Identity();
    }
    else if (almost_equal(lat.c(), -raw_reduced[2]))
    {
        T = Eigen::Matrix3d::Identity() * S;
    }

    else if (almost_equal(lat.c(), raw_reduced[1]))
    {
        T = P;
    }
    else if (almost_equal(lat.c(), -raw_reduced[1]))
    {
        T = P * S;
    }

    else if (almost_equal(lat.c(), raw_reduced[0]))
    {
        T = P * P;
    }
    else if (almost_equal(lat.c(), -raw_reduced[0]))
    {
        T = P * P * S;
    }

    else
    {
        auto [V, E] = approximate_integer_transformation(lat.column_vector_matrix(), raw_reduced.lat_column_mat());
        throw std::runtime_error("Could not permute lattice vectors to recover C after making reduced cell");
    }

    return xtal::Lattice::from_column_vector_matrix(raw_reduced.lat_column_mat() * T);
}

MoireGenerator::MoireGenerator(const xtal::Lattice& input_lat, double degrees, long max_lattice_sites)
    : moire(input_lat, degrees),
      transformation_matrix_to_super_aligned_moire(Eigen::Matrix3l::Identity()),
      transformation_matrix_to_super_rotated_moire(Eigen::Matrix3l::Identity()),
      aligned_moire_approximant(moire.aligned_moire_lattice, moire.aligned_lattice, moire.rotated_lattice),
      rotated_moire_approximant(moire.rotated_moire_lattice, moire.aligned_lattice, moire.rotated_lattice)
{
    // Maps to make a for loop over aligned and rotated memebers
    std::unordered_map<LATTICE, const xtal::Lattice*> moire_units;
    moire_units.emplace(LATTICE::ALIGNED, &moire.aligned_moire_lattice);
    moire_units.emplace(LATTICE::ROTATED, &moire.rotated_moire_lattice);

    std::unordered_map<LATTICE, MoireApproximant*> best_approximants;
    best_approximants.emplace(LATTICE::ALIGNED, &this->aligned_moire_approximant);
    best_approximants.emplace(LATTICE::ROTATED, &this->rotated_moire_approximant);

    Eigen::Matrix3l tmp;
    std::unordered_map<LATTICE, Eigen::Matrix3l*> trans_mats_to_super_moire;
    trans_mats_to_super_moire.emplace(LATTICE::ALIGNED, &this->transformation_matrix_to_super_aligned_moire);
    trans_mats_to_super_moire.emplace(LATTICE::ROTATED, &this->transformation_matrix_to_super_rotated_moire);

    // Repeat the steps for things generated from the aligned and rotated zones
    for (ZONE lat : {ZONE::ALIGNED, ZONE::ROTATED})
    {
        const auto& moire_unit = *moire_units[lat];
        const auto& best_approximant = *best_approximants[lat];

        long min_lattice_sites =
            best_approximant.approximate_moire_integer_transformations.at(LATTICE::ALIGNED).determinant() +
            best_approximant.approximate_moire_integer_transformations.at(LATTICE::ROTATED).determinant();

        long max_moire_scel_size = max_lattice_sites / min_lattice_sites;

        // If the max allowed number of Moirons is less than two, don't bother trying to find supercells
        if (max_moire_scel_size < 2)
        {
            continue;
        }

        CASM::xtal::ScelEnumProps super_moire_props(2, max_moire_scel_size + 1, "ab");
        CASM::xtal::SuperlatticeEnumerator super_moire_enumerator(
            moire_unit.__get(), xtal::make_point_group(moire_unit, 1e-5), super_moire_props);

        const auto& aligned_unit = this->moire.aligned_lattice;
        const auto& rotated_unit = this->moire.rotated_lattice;

        double best_error = error_metric(best_approximant.approximate_moire_lattice, aligned_unit, rotated_unit);

        // For each possible Moire supercell, see if it coincides better, if so, keep it.
        for (auto super_moire_it = super_moire_enumerator.begin(); super_moire_it != super_moire_enumerator.end();
             ++super_moire_it)
        {
            xtal::Lattice super_moire = make_reduced_cell(*super_moire_it);
            // the C vectors should't be changing
            assert(CASM::almost_equal(super_moire.c(), super_moire_it->operator[](2).eval()));

            double current_error = error_metric(super_moire, aligned_unit, rotated_unit);

            auto [moire_to_super_trans_mat, residual] = approximate_integer_transformation(moire_unit, super_moire);
            assert(almost_zero(residual));

            if (current_error < best_error-1e-5)
            {
                best_error = current_error;
                *best_approximants[lat] = MoireApproximant(super_moire, aligned_unit, rotated_unit);
                *trans_mats_to_super_moire[lat] = moire_to_super_trans_mat;
            }
        }
    }
}

double
MoireGenerator::error_metric(const xtal::Lattice& moire, const xtal::Lattice& aligned, const xtal::Lattice& rotated)
{
    MoireApproximant approx(moire,aligned,rotated);
    double error = 0.0;

    for (auto LAT : {LATTICE::ALIGNED, LATTICE::ROTATED})
    {
        for (int i : {0, 1})
        {
            if(!almost_zero(approx.approximate_moire_integer_transformations[LAT]))
            {
                error += approx.approximate_moire_integer_transformation_errors[LAT].col(i).norm();
            }
        }

        //This works really poorly
        /* const auto& F=approx.approximation_deformations[LAT]; */
        /* DeformationReport report(F); */
        /* error+=report.deviatoric_strain; */
    }
    return error;
}

MoireStructureGenerator::MoireStructureGenerator(const Structure& slab_unit, double degrees, long max_lattice_sites)
    : MoireGenerator(slab_unit.lattice(), degrees, max_lattice_sites), slab_unit(slab_unit)
{
}

MoireStructureGenerator::Structure MoireStructureGenerator::layer(ZONE brillouin, LATTICE lat) const
{
    const auto report=this->generate(brillouin,lat);

    const auto& approx_lat=report.approximate_tiling_unit;
    const auto& T= report.tiling_unit_supercell_matrix;
    /* const auto& approx_lat = this->approximate_lattice(brillouin, lat); */
    /* const auto& T = this->approximate_moire_integer_transformation(brillouin, lat); */

    Structure approx_unit = slab_unit;
    approx_unit.set_lattice(approx_lat, xtal::FRAC);

    auto layer = xtal::make_superstructure(approx_unit, T.cast<int>());
    return layer;
}
} // namespace mush
} // namespace casmutils
