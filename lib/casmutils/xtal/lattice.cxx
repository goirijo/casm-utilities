#include <casm/crystallography/Niggli.hh>
#include <casm/misc/CASM_Eigen_math.hh>
#include <casmutils/xtal/lattice.hpp>

namespace
{
using namespace casmutils;
/// Ensure that both lattices have parallel vectors axb
bool ab_plane_conserved(const xtal::Lattice& lhs, const xtal::Lattice& rhs)
{
    Eigen::Vector3d lhs_norm = lhs[0].cross(lhs[1]).normalized();
    Eigen::Vector3d rhs_norm = rhs[0].cross(rhs[1]).normalized();

    double norm_dot = lhs_norm.dot(rhs_norm);
    if (CASM::almost_equal(norm_dot, 1.0) || CASM::almost_equal(norm_dot, -1.0))
    {
        return true;
    }

    return false;
}
} // namespace

namespace casmutils
{
namespace xtal
{
Lattice::Lattice(const CASM::xtal::Lattice& init_lat) : casm_lattice(init_lat) {}
Lattice::Lattice(const Eigen::Matrix3d& column_lat_mat) : casm_lattice(column_lat_mat) {}

Lattice::Lattice(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
    : Lattice(Lattice::stack_column_vectors(a, b, c))
{
}

Eigen::Matrix3d
Lattice::stack_column_vectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    Eigen::Matrix3d column_matrix;
    column_matrix.col(0) = a;
    column_matrix.col(1) = b;
    column_matrix.col(2) = c;
    return column_matrix;
}

LatticeEquals_f::LatticeEquals_f(const Lattice& ref_lat, double tol) : ref_lat(ref_lat), tol(tol) {}
bool LatticeEquals_f::operator()(const Lattice& other)
{
    return ref_lat.column_vector_matrix().isApprox(other.column_vector_matrix(), tol);
}

void make_niggli(Lattice* lattice_ptr)
{
    *lattice_ptr = CASM::xtal::niggli(CASM::xtal::Lattice(lattice_ptr->column_vector_matrix()), CASM::TOL);
    return;
}
Lattice make_niggli(const Lattice& non_niggli_lattice)
{
    Lattice niggli_lattice = non_niggli_lattice;
    make_niggli(&niggli_lattice);
    return niggli_lattice;
}

Lattice make_reciprocal(const Lattice& real_lattice)
{
    auto casm_reciprocal = real_lattice.__get().reciprocal();
    // TODO: Move constructor for lattice?
    return Lattice(casm_reciprocal);
}

std::pair<Eigen::Matrix3d, Eigen::Matrix3d> polar_decomposition(Eigen::Matrix3d const& F)
{
    Eigen::Matrix3d U = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(F * F.transpose()).operatorSqrt();
    Eigen::Matrix3d R = U.inverse() * F;
    return std::make_pair(R, U);
}

Lattice slice_along_plane(const Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes)
{
    // TODO: Fails if you give it stupid values like (0,-2,0);
    auto orthoscore = [=](const Lattice& l) {
        return std::abs(l.a().normalized().dot(l.b().normalized()));
    };

    // The 0 means "get the smallest cell possible", and I wish it was the default
    Lattice shift_lattice = unit_lattice.__get().lattice_in_plane(miller_indexes, 0);

    Lattice best_lattice = shift_lattice;
    double last_score = orthoscore(shift_lattice);
    for (const Eigen::Matrix3i& mat : CASM::unimodular_matrices())
    {
        // discard any transformation on the a or b vector that isn't a linear combination of a and b
        // discard any transformation that modifies c
        if (mat(0, 2) != 0 || mat(1, 2) != 0 || mat(2, 0) != 0 || mat(2, 1) != 0 || mat(2, 2) != 1)
        {
            continue;
        }

        // discard trasnformations that switch the handedness of the lattice
        if (mat.determinant() * unit_lattice.column_vector_matrix().determinant() < 0)
        {
            continue;
        }

        /* auto candidate_lat_mat = shift_lattice.column_vector_matrix() * mat.cast<double>(); */
        Lattice new_lattice = make_superlattice(shift_lattice, mat);
        double new_score = orthoscore(new_lattice);
        if (new_score < last_score)
        {
            best_lattice = new_lattice;
            last_score = new_score;
        }
    }

    assert(::ab_plane_conserved(best_lattice, shift_lattice));
    return best_lattice;
}

Lattice make_right_handed(const Lattice& left_handed_lattice)
{
    auto casm_lattice = left_handed_lattice.__get();
    casm_lattice.make_right_handed();
    return Lattice(casm_lattice);
}

Lattice make_superlattice(const Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat)
{
    return Lattice(CASM::xtal::make_superlattice(tiling_unit.__get(), col_transf_mat));
}

} // namespace xtal
} // namespace casmutils
