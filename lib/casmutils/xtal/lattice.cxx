#include <casm/crystallography/Niggli.hh>
#include <casmutils/xtal/lattice.hpp>
#include <cmath>

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

Lattice Lattice::from_lattice_parameters(
    const double& a, const double& b, const double& c, const double& alpha, const double& beta, const double& gamma)
{
    // Convert angles to radians
    double alpha_rad = (alpha * M_PI) / 180.;
    double beta_rad = (beta * M_PI) / 180.;
    double gamma_rad = (gamma * M_PI) / 180.;

    // Compute sins and cosines
    double cos_alpha = std::cos(alpha_rad);
    double cos_beta = std::cos(beta_rad);
    double cos_gamma = std::cos(gamma_rad);
    double sin_alpha = std::sin(alpha_rad);
    double sin_beta = std::sin(beta_rad);
    double sin_gamma = std::sin(gamma_rad);

    Eigen::Vector3d a_vec{a * sin_beta, 0, a * cos_beta};
    Eigen::Vector3d c_vec{0, 0, c};

    double b_x = (cos_gamma - cos_alpha * cos_beta) / sin_beta;
    double b_y = std::sqrt((sin_alpha * sin_alpha) - (b_x * b_x));
    Eigen::Vector3d b_vec{b * b_x, b * b_y, b * cos_alpha};

    return Lattice(a_vec, b_vec, c_vec);
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

} // namespace xtal
} // namespace casmutils
