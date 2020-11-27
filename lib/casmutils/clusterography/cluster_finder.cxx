#include <casm/clusterography/ClusterInvariants.hh>
#include <casm/clusterography/ClusterSpecs.hh>
#include <casm/clusterography/ClusterSymCompareDecl.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/symmetry/Orbit.hh>
#include <casm/symmetry/OrbitDecl.hh>
#include <casmutils/clusterography/cluster_finder.hpp>
#include <functional>
#include <memory>

bool return_true(CASM::xtal::Site site) { return true; }

void make_periodic_orbits(std::vector<double> max_length, const casmutils::xtal::Structure structure)
{
    CASM::xtal::BasicStructure basic_structure = structure.__get<CASM::xtal::BasicStructure>();
    const CASM::Structure struc(basic_structure);
    auto shared_prim = std::make_shared<const CASM::Structure>(struc);
    auto sym_group = shared_prim->factor_group();
    std::function<bool(CASM::xtal::Site)> site_filter = return_true;

    CASM::PeriodicMaxLengthClusterSpecs clus_specs(shared_prim, sym_group, site_filter, max_length);
    CASM::ClusterSpecs::PeriodicOrbitVec orbits = clus_specs.make_periodic_orbits(std::cout);
    std::cout << orbits.size() << std::endl;
    std::cout << "---------" << std::endl;
    std::cout << orbits[16][1].size() << std::endl;
}
