#ifndef CLUSTER_FINDER_HH
#define CLUSTER_FINDER_HH

#include <casm/crystallography/Structure.hh>
#include <casmutils/sym/cartesian.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure.hpp>

/// File contents:
/// This file provides structs & functions related to finding clusters and orbits for a structure:
/// To be able to compile this you need to compile casmutils with --enable-fullcasm option since this module
/// utilizes symmetry & clustergraphy modules of casm
/// In namespace - casutils::cluster,
/// structs available: Cluster, Orbit - Look for the definitions below to use them
/// functions available: make_periodic_clusters, make_local_clusters - Look for the definitions below to use them

namespace casmutils
{
namespace cluster
{

/// A cluster is defined as a group of sites
/// TODO: Make it available without --enable-fullcasm
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
/// TODO: Make it available without --enable-fullcasm
/// TODO: It still has only basic functionality but should be improved
/// TODO: Get the symmetry information - Generating Group from prototype
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

    // Returns the prototype of the orbit - Equivalent to 0th element in the orbit
    Cluster prototpye() const { return this->m_clusters[0]; }

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

/// Default site_filter function
bool default_site_filter(const xtal::Site& site);

/// Returns periodic orbits of a given structure - PeriodicOrbits suggests that the clusters in any orbit branch
/// need not be in the first unit cell and it entirely depends on the max_length parameter which is explained below
/// @params - max_length[b] : Includes clusters in an orbit branch "b" only if the distance between any two sites
/// is less than max_length[b]. For example, Orbit branch "2" indicates all the 2 point clusters, orbit branch "3"
/// indicates all the 3 point clusters.
/// The first and the second elements in the max_length vector will not be considered (since they belong to
/// orbit branch 0 & 1 and max_length does not make sense for null point & single point clusters)
/// A sample max_length can be {2, 2, 3, 4} - you will get upto 3 point clusters with max_length of 2 point
/// clusters is 3 Angstroms & 3 point clusters is 4 Angstroms.
/// A caveat - The algorithm is designed in such a way that it uses the previous orbit branch to figure out the
/// new orbit branch. For example, in the above mentioned case of max_length, the 3 point cluster orbits will be
/// generated using the 2 point cluster orbits. Hence, if you want all the 3 point clusters with max_length 4,
/// you are required to have max_length for 2 point clusters to be at least 4.
/// @params - structure - casmutils::xtal::Structure for which you are finding clusters
/// @params - site_filter is a function which when provided casmutils::xtal::Site as an argument returns a boolean
/// true or false. A default site filter function is provided already where all the sites by default are included.
std::vector<Orbit> make_periodic_orbits(const std::vector<double>& max_length,
                                        const casmutils::xtal::Structure& structure,
                                        const std::function<bool(const xtal::Site)>& site_filter = default_site_filter);

/// Returns local orbits of a given structure - LocalOrbits suggests that the clusters in any orbit branch will
/// always be within a cut off radius of a given phemomenal cluster
/// @params - phenomenal_cluster - It is a casmutils::Cluster object around which you want your local clusters
/// @params - include_phenomenal_cluster - A boolean flag which suggests if the phenomenal cluster should be
/// included in your output. If true, the phenomenal cluster will be added as a seperate orbit. It won't be
/// grouped together with equivalent clusters of the phenomenal cluster. They will belong to a different orbit
/// @params - max_length[b] : Includes clusters in an orbit branch "b" only if the distance between any two sites
/// is less than max_length[b]. For example, Orbit branch "2" indicates all the 2 point clusters, orbit branch "3"
/// indicates all the 3 point clusters.
/// The first and the second elements in the max_length vector will not be considered (since they belong to
/// orbit branch 0 & 1 and max_length does not make sense for null point & single point clusters)
/// A sample max_length can be {2, 2, 3, 4} - you will get upto 3 point clusters with max_length of 2 point
/// clusters is 3 Angstroms & 3 point clusters is 4 Angstroms.
/// @params - cut_off_radius[b] - Includes clusters in an orbit branch "b" if all the sites are within the
/// provided cut_off_radius of any of the sites in the phenomenal cluster. The first element of the vector is
/// always ignored as cut_off_radius doesn't make sense for a null point cluster
/// A caveat - The algorithm is designed in such a way that it uses the previous orbit branch to figure out the
/// new orbit branch. For example, in the above mentioned case of max_length, the 3 point cluster orbits will be
/// generated using the 2 point cluster orbits. Hence, if you want all the 3 point clusters with max_length 4,
/// you are required to have max_length for 2 point clusters to be at least 4. This also applies for
/// cut_off_radius.
/// @params - structure - casmutils::xtal::Structure for which you are finding clusters
/// @params - site_filter is a function which when provided casmutils::xtal::Site as an argument returns a boolean
/// true or false. A default site filter function is provided already where all the sites by default are included.
std::vector<Orbit> make_local_orbits(const Cluster& phenomenal_cluster,
                                     bool include_phenomenal_cluster,
                                     const std::vector<double>& max_length,
                                     const std::vector<double>& cut_off_radius,
                                     const casmutils::xtal::Structure& structure,
                                     const std::function<bool(const xtal::Site)>& site_filter = default_site_filter);

} // namespace cluster
} // namespace casmutils

namespace extend
{
typedef casmutils::cluster::Cluster Cluster;
typedef casmutils::cluster::Orbit Orbit;

template <typename CasmOrbitType>
Orbit casm_orbit_to_casmutils_orbit(const CasmOrbitType& casm_orbit, const CASM::Structure& casm_struc)
{
    std::vector<Cluster> clusters;
    for (const auto& casm_cluster : casm_orbit)
    {
        Cluster casm_utils_cluster;
        for (const auto& casm_element : casm_cluster.elements())
        {
            auto casm_site = casm_element.site(casm_struc);
            casm_utils_cluster.add_site(extend::casm_site_to_casmutils_site(casm_site));
        }

        clusters.push_back(casm_utils_cluster);
    }
    return Orbit(clusters);
}
} // namespace extend
#endif /* ifndef CLUSTER_FINDER_H */
