#include <casm/clusterography/ClusterInvariants.hh>
#include <casm/clusterography/ClusterSpecs.hh>
#include <casm/clusterography/ClusterSymCompareDecl.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/symmetry/Orbit.hh>
#include <casm/symmetry/OrbitDecl.hh>
#include <casmutils/clusterography/cluster_finder.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <functional>
#include <memory>

namespace casmutils
{
namespace clusterography
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

bool return_true(CASM::xtal::Site site) { return true; }

std::vector<Orbit> make_periodic_orbits(const std::vector<double>& max_length,
                                        const casmutils::xtal::Structure& structure)
{
    CASM::xtal::BasicStructure basic_structure = structure.__get<CASM::xtal::BasicStructure>();
    const CASM::Structure struc(basic_structure);
    auto shared_prim = std::make_shared<const CASM::Structure>(struc);
    auto sym_group = shared_prim->factor_group();
    std::function<bool(CASM::xtal::Site)> site_filter = return_true;

    CASM::PeriodicMaxLengthClusterSpecs clus_specs(shared_prim, sym_group, site_filter, max_length);
    CASM::ClusterSpecs::PeriodicOrbitVec casm_orbits = clus_specs.make_periodic_orbits(std::cout);

    std::vector<Orbit> orbits;

    for (const auto& casm_orbit : casm_orbits)
    {
        std::vector<Cluster> clusters;

        for (const auto& casm_cluster : casm_orbit)
        {
            Cluster casm_utils_cluster;
            for (const auto& casm_element : casm_cluster.elements())
            {
                auto casm_site = casm_element.site(*shared_prim);
                Eigen::Vector3d cart_coords = casm_site.cart();
                std::string label;
                for (const auto& m : casm_site.occupant_dof())
                {
                    label += m.name();
                }

                xtal::Site casm_utils_site(cart_coords, label);
                casm_utils_cluster.add_site(casm_utils_site);
            }

            clusters.push_back(casm_utils_cluster);
        }

        orbits.push_back(Orbit(clusters));
    }

    return orbits;
}

} // namespace clusterography
} // namespace casmutils
