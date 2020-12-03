#include <casm/clusterography/ClusterInvariants.hh>
#include <casm/clusterography/ClusterSpecs.hh>
#include <casm/clusterography/ClusterSymCompare.hh>
#include <casm/clusterography/ClusterSymCompareDecl.hh>
#include <casm/clusterography/IntegralCluster.hh>
#include <casm/crystallography/UnitCellCoord.hh>
#include <casm/symmetry/InvariantSubgroup.hh>
#include <casm/symmetry/Orbit.hh>
#include <casm/symmetry/OrbitDecl.hh>
#include <casmutils/clusterography/cluster_finder.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <cassert>
#include <functional>
#include <memory>

namespace casmutils
{
namespace cluster
{

xtal::Site Cluster::operator[](const std::size_t& i) const
{
    if (i >= this->m_sites.size())
    {
        throw std::out_of_range("Site out of range");
    }
    return this->m_sites[i];
}

Eigen::Vector3d Cluster::geometric_center_of_mass() const
{
    Eigen::Vector3d center_of_mass;

    for (const auto& site : this->m_sites)
    {
        center_of_mass += site.cart();
    }

    return center_of_mass / this->size();
}

Cluster Orbit::operator[](const std::size_t& i) const
{
    if (i >= this->m_clusters.size())
    {
        throw std::out_of_range("Cluster out of range");
    }
    return this->m_clusters[i];
}

Cluster operator*(const sym::CartOp& sym_op, const Cluster& cluster)
{
    std::vector<xtal::Site> new_sites;

    for (const auto& site : cluster)
    {
        new_sites.push_back(sym_op * site);
    }

    return Cluster(new_sites);
}

bool default_site_filter(const xtal::Site& site) { return true; }

std::vector<Orbit> make_periodic_orbits(const std::vector<double>& max_length,
                                        const casmutils::xtal::Structure& structure,
                                        const std::function<bool(const xtal::Site)>& site_filter)
{
    // Get a CASM::Structure ptr
    CASM::xtal::BasicStructure basic_structure = structure.__get<CASM::xtal::BasicStructure>();
    const CASM::Structure casm_struc(basic_structure);
    auto shared_prim = std::make_shared<const CASM::Structure>(casm_struc);

    // Generating SymGroup for clusters
    auto sym_group = shared_prim->factor_group();

    // Make casm SiteFilterFunction from casm utils site_filter
    std::function<bool(CASM::xtal::Site)> casm_site_filter = [&](const CASM::xtal::Site& casm_site) {
        return site_filter(extend::casm_site_to_casmutils_site(casm_site));
    };

    // Make casm_orbits
    CASM::PeriodicMaxLengthClusterSpecs clus_specs(shared_prim, sym_group, casm_site_filter, max_length);
    CASM::ClusterSpecs::PeriodicOrbitVec casm_orbits = clus_specs.make_periodic_orbits(std::cout);

    // Convert to casmutils orbits
    std::vector<Orbit> orbits;
    for (const auto& casm_orbit : casm_orbits)
    {
        orbits.push_back(extend::casm_orbit_to_casmutils_orbit(casm_orbit, *shared_prim));
    }

    return orbits;
}

std::vector<Orbit> make_local_orbits(const Cluster& phenomenal_cluster,
                                     bool include_phenomenal_cluster,
                                     const std::vector<double>& max_length,
                                     const std::vector<double>& cut_off_radius,
                                     const casmutils::xtal::Structure& structure,
                                     const std::function<bool(const xtal::Site)>& site_filter)
{
    // Assert the size of max_length & cut_off radius to be the same. Otherwise it is undefined behavior
    assert((max_length.size() == cut_off_radius.size()) &&
           "Orbit branches to be considered should be the same in max_length & cut_off_radius");

    // Get a CASM::Structure ptr
    CASM::xtal::BasicStructure basic_structure = structure.__get<CASM::xtal::BasicStructure>();
    const CASM::Structure casm_struc(basic_structure);
    auto shared_prim = std::make_shared<const CASM::Structure>(casm_struc);

    // Get the SiteFilterFunction for casm from casmutils site_filter
    std::function<bool(CASM::xtal::Site)> casm_site_filter = [&](const CASM::xtal::Site& casm_site) {
        return site_filter(extend::casm_site_to_casmutils_site(casm_site));
    };

    // Convert the phenomenal cluster to CASM::IntegralCluster
    std::vector<CASM::xtal::UnitCellCoord> unit_cell_coords;
    for (const auto& site : phenomenal_cluster)
    {
        CASM::xtal::UnitCellCoord unit_cell_coord = CASM::xtal::UnitCellCoord::from_coordinate(
            shared_prim->structure(),
            CASM::xtal::Coordinate(site.cart(), shared_prim->structure().lattice(), CASM::CART),
            CASM::TOL);
        unit_cell_coords.push_back(unit_cell_coord);
    }
    CASM::IntegralCluster casm_phenomenal_cluster(*shared_prim, unit_cell_coords.begin(), unit_cell_coords.end());

    // Get the invariant SymGroup of the phenomenal cluster
    CASM::LocalSymCompare<CASM::IntegralCluster> local_sym_compare(shared_prim, shared_prim->lattice().tol());
    auto sym_group =
        CASM::make_invariant_subgroup(casm_phenomenal_cluster, shared_prim->factor_group(), local_sym_compare);

    // Make casm_orbits
    CASM::LocalMaxLengthClusterSpecs cluster_specs(shared_prim,
                                                   sym_group,
                                                   casm_phenomenal_cluster,
                                                   casm_site_filter,
                                                   max_length,
                                                   cut_off_radius,
                                                   include_phenomenal_cluster);

    CASM::ClusterSpecs::LocalOrbitVec casm_orbits = cluster_specs.make_local_orbits(std::cout);

    // Convert casm_orbits to casmutils orbits
    std::vector<Orbit> orbits;
    for (const auto& casm_orbit : casm_orbits)
    {
        orbits.push_back(extend::casm_orbit_to_casmutils_orbit(casm_orbit, *shared_prim));
    }

    return orbits;
}

} // namespace cluster
} // namespace casmutils
