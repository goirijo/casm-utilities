#ifndef CLUSTER_FINDER_HH
#define CLUSTER_FINDER_HH

#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>

namespace casmutils
{
namespace clusterography
{

/// A cluster is defined as a group of sites
struct Cluster
{
public:
    // Constructors - Cluster can be constructed with a vector of sites
    Cluster() {} // Default constructor - Implies an empty cluster
    Cluster(const std::vector<xtal::Site>& sites) : m_sites(sites) {}

    // Returns the sites in the cluster
    std::vector<xtal::Site> get_sites() const { return this->m_sites; }

    // Returns the number of sites in the cluster
    int size() const { return this->m_sites.size(); }

    // Const Iterators - Cluster can be used in range based for loop
    std::vector<xtal::Site>::const_iterator begin() const { return this->m_sites.begin(); }
    std::vector<xtal::Site>::const_iterator end() const { return this->m_sites.end(); }

    // Access element i in the cluster
    xtal::Site operator[](const std::size_t& i) const;

    // Calculates the geometric center of mass and returns the value in cartesian coordinates
    Eigen::Vector3d geometric_center_of_mass() const;

    // Adds a site to the cluster
    void add_site(const xtal::Site& site) { this->m_sites.push_back(site); }
    // TODO: Any other useful methods for a cluster?

private:
    std::vector<xtal::Site> m_sites;
};

/// An orbit is defined as a group of symmetrically equivalent clusters
/// TODO: It still has only basic functionality but should be improved
struct Orbit
{
public:
    // Constructor - Currently Orbit can be constructed with a vector of clusters
    // TODO: Orbit by definition contains symmetrically equivalent clusters. If you have a constructor
    // with just a vector of clusters it defeats the purpose. Need to do something about it
    // TODO: Need to have a cluster compare function which can take care of the above mentioned issue?
    // TODO: Can you have an incomplete orbit as well?
    // TODO: Can have a constructor which can generate the orbit directly from a prototype and symmetry info
    Orbit(const std::vector<Cluster>& clusters) : m_clusters(clusters) {}

    // Returns all the clusters in the orbit
    std::vector<Cluster> get_clusters() const { return this->m_clusters; }

    // Returns the number of clusters in the orbit
    int size() const { return this->m_clusters.size(); }

    // Const iterators - Orbit can be used in range based for loop
    std::vector<Cluster>::const_iterator begin() const { return this->m_clusters.begin(); }
    std::vector<Cluster>::const_iterator end() const { return this->m_clusters.end(); }

    // Access cluster i in the orbit
    Cluster operator[](const std::size_t& i) const;

    // TODO: Any useful methods?

private:
    std::vector<Cluster> m_clusters;
};

/// Apply a given symmetry operation to a cluster
Cluster operator*(const sym::CartOp& sym_op, const Cluster& cluster);

/// Returns periodic orbits of a given structure
/// TODO: Explain args further and include SiteFilterFunction as an arg
std::vector<Orbit> make_periodic_orbits(const std::vector<double>& max_length,
                                        const casmutils::xtal::Structure& structure);

// TODO: Write a function which does local orbits

} // namespace clusterography
} // namespace casmutils

#endif /* ifndef CLUSTER_FINDER_H */
