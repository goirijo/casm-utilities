#ifndef ROCKSALTTOGGLER_HH

#include <string>
#include <unordered_map>
#include <utility>

#include <casmutils/definitions.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>

namespace enumeration
{
class RockSaltOctahedraToggler
{
public:
    typedef Eigen::Vector3d Coordinate;
    typedef casmutils::xtal::Structure Structure;
    typedef casmutils::xtal::Lattice Lattice;
    typedef int index;

    /// Construct with a transformation matrix relative to the primitive structure.
    /// The species names refers to the species of the cation and anions.
    /// The central specie specifies wheter the cation or anion is at the origin, and will
    /// be the central atom of the octahedra.
    static RockSaltOctahedraToggler relative_to_primitive(const Eigen::Matrix3i trans_mat,
                                                          std::pair<std::string, std::string> species_names,
                                                          bool second_is_central,
                                                          double init_nn_distance);

    // TODO
    /// Construct with a transformation matrix relative to the conventional structure
    /// The species names refers to the species of the cation and anions.
    /// The central specie specifies wheter the cation or anion is at the origin, and will
    /// be the central atom of the octahedra.
    static RockSaltOctahedraToggler relative_to_conventional(const Eigen::Matrix3i trans_mat,
                                                             std::pair<std::string, std::string> species_names,
                                                             bool second_is_central,
                                                             double init_nn_distance);

    /// Sets the central coordinate (center of octahedron) ON
    /// along with its surrounding oxygen (or other specified anion)
    void activate(const Coordinate& central_coord);
    void activate(index central_ix);
    /// Sets all octahera ON
    void activate_all();

    /// Sets the central coordinate (center of octahedron ) OFF
    /// along with its surrounding oxygen (or other specified anion),
    /// if they are not part of another octahedron
    void deactivate(const Coordinate& central_coord);
    void deactivate(index central_ix);
    /// Sets all octahera OFF
    void deactivate_all();

    /// Calls activate/deactivate to reverse whether the octahedron is there or not
    void toggle(const Coordinate& central_coord);
    void toggle(index central_ix);
    /// Calls activate/deactivate to reverse whether the octahedron is there or not on every octahedron
    void toggle_all();

    /// Print the current state of the structure to a stream
    void print(std::ostream& out_stream) const;

    /// Get list of available coordinates for centers of possible octahedra (cation sites)
    std::vector<std::pair<index, Coordinate>> all_octahedron_center_coordinates() const;

    /// Get a representation of the rocksalt structure in its current state
    Structure structure() const;

    /// Return the structure for the primitive rocksalt
    static Structure primitive_structure(std::pair<std::string, std::string> species_names, double init_nn_distance);

    /// Return the structure for the conventional cell of  rocksalt
    static Structure conventional_structure(std::pair<std::string, std::string> species_names, double init_nn_distance);

    /// Defines the nearest neighbor distance (distance between central ion and vertex)
    double nearest_neighbor_distance() const;

private:
    /// Defines the nearest neighbor distance (distance between central ion and vertex)
    double nn_distance;

    /// The rocksalt structure we're working with, specified at construction by transformation matrix
    /// and species
    mutable Structure rocksalt_struc;

    /// This is the species that makes the corners of the octahedra
    const std::string vertex_ion_name;

    /// This is the species at the center of the octahedra
    const std::string central_ion_name;

    /// Distances to each of the nearest neighbor atoms from the center of the octahedron
    const std::array<Coordinate, 6> nearest_neighbor_deltas;

    /// Keeps track of which central_ions are turned on
    std::unordered_map<index, bool> central_ion_is_on;

    /// Counts by how many central_ions an vertex_ion is being held by
    /// Useful to know when the vertex_ion can be "released" and switched off (when count becomes zero)
    std::unordered_map<index, int> leashed_vertex_ions;

    RockSaltOctahedraToggler(Structure&& init_struc,
                             std::string init_central_name,
                             std::string init_vertex_name,
                             std::array<Coordinate, 6>&& init_nn_deltas,
                             std::unordered_map<index, bool>&& init_central_is_on,
                             std::unordered_map<index, int>&& init_leashes,
                             double init_nn_distance);

    /// Converts the given coordinate to the corresponding index within the Structure
    index coordinate_to_index(Coordinate coordinate) const;
    Coordinate index_to_coordinate(index coordinate_index) const;

    /// Retruns true if the given index exists as a possible central_ion site in the rocksalt structure
    bool is_central_ion_index(index ix) const;

    /// Retruns true if the given index exists as a possible vertex_ion site in the rocksalt structure
    bool is_vertex_ion_index(index ix) const;

    ///
    void turn_central_ion_on(index ix);
    void turn_central_ion_off(index ix);

    /// Returns coordinates of the nearest neighbors for the given site index
    std::array<Coordinate, 6> nearest_neighbor_site_coordinates(const Coordinate& central_coord) const;
    std::array<index, 6> nearest_neighbor_site_indexes(index central_coord_index) const;

    /// For each of the given vertex_ion coordinates, increment their leashed count by 1
    void increment_leashes(const std::array<Coordinate, 6>& neighboring_vertex_ion_coords);

    /// For each of the given vertex_ion coordinates, decrease their leashed count by 1
    void reduce_leashes(const std::array<Coordinate, 6>& neighboring_vertex_ion_coords);

    /// Apply the current state of the octahedra to the structure
    void commit() const;
    void commit_central_ions() const;
    void commit_vertex_ions() const;

    // helper functions for constructor
    static double primitive_lattice_scale_factor();
    static std::array<Coordinate, 6> initialized_nearest_neighbor_deltas(double init_nn_distance);
    static std::unordered_map<index, bool> initialized_central_ion_is_on(const Structure& init_struc,
                                                                         std::string central_name);
    static std::unordered_map<index, int> initialized_leashed_vertex_ions(const Structure& init_struc,
                                                                          std::string vertex_name);
};
} // namespace enumeration

#endif
