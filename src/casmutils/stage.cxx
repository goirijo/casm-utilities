#include "casmutils/stage.hpp"
#include "casm/crystallography/Structure.hh"
#include "casmutils/exceptions.hpp"

namespace
{
// Please never use this outisde of the RockSalt context
void set_site_occupant(CASM::Site* mutating_site, std::string new_occ)
{
    if (new_occ == "Va")
    {
        mutating_site->set_occ_value(0);
        assert(mutating_site->allowed_occupants()[0] == "Va");
    }

    else
    {
        mutating_site->set_occ_value(1);
    }
}
} // namespace

// TODO: CAPS??
namespace Extend
{
CASM::Site multi_atomic_site(const CASM::Coordinate& coord, const std::vector<std::string>& allowed_species)
{
    CASM::Array<CASM::Molecule> allowed_molecules;
    for (auto specie : allowed_species)
    {
        //make_atom returns a molecule lol
        allowed_molecules.push_back(make_atom(specie, coord.home()));
    }

    CASM::Site site(coord, "");
    site.set_site_occupant(allowed_molecules);

    return site;
}
} // namespace Extend

namespace SpecializedEnumeration
{

void RockSaltOctahedraToggler::activate(index central_coord_index)
{
    if (!this->is_central_ion_index(central_coord_index))
    {
        throw UtilExcept::IncompatibleCoordinate();
    }

    if (central_ion_is_on.at(central_coord_index))
    {
        return;
    }

    this->turn_central_ion_on(central_coord_index);

    auto nearest_vertex_ions = this->nearest_neighbor_site_coordinates(this->index_to_coordinate(central_coord_index));
    this->increment_leashes(nearest_vertex_ions);

    return;
}

void RockSaltOctahedraToggler::activate(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);
    this->activate(central_coord_index);
    return;
}

void RockSaltOctahedraToggler::deactivate(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);
    if (!this->is_central_ion_index(central_coord_index))
    {
        throw UtilExcept::IncompatibleCoordinate();
    }

    if (!central_ion_is_on.at(central_coord_index))
    {
        return;
    }

    this->turn_central_ion_off(central_coord_index);

    auto nearest_vertex_ions = this->nearest_neighbor_site_coordinates(central_coord);
    this->reduce_leashes(nearest_vertex_ions);

    return;
}

bool RockSaltOctahedraToggler::is_central_ion_index(index ix) const
{
    return central_ion_is_on.find(ix) != central_ion_is_on.end();
}
bool RockSaltOctahedraToggler::is_vertex_ion_index(index ix) const
{
    return leashed_vertex_ions.find(ix) != leashed_vertex_ions.end();
}

void RockSaltOctahedraToggler::turn_central_ion_on(index ix) { central_ion_is_on[ix] = true; }
void RockSaltOctahedraToggler::turn_central_ion_off(index ix) { central_ion_is_on[ix] = false; }

std::array<RockSaltOctahedraToggler::Coordinate, 6>
RockSaltOctahedraToggler::nearest_neighbor_site_coordinates(const Coordinate& central_coord) const
{
    std::array<Coordinate, 6> nearest_coords{central_coord, central_coord, central_coord,
                                             central_coord, central_coord, central_coord};
    for (int i = 0; i < 6; ++i)
    {
        nearest_coords[i] += nearest_neighbor_deltas[i];
    }

    return nearest_coords;
}

std::array<RockSaltOctahedraToggler::index, 6>
RockSaltOctahedraToggler::nearest_neighbor_site_indexes(index central_coord_index) const
{
    auto nearest_coordinates = this->nearest_neighbor_site_coordinates(this->index_to_coordinate(central_coord_index));
    std::array<index, 6> nearest_coord_indexes;

    for (int i = 0; i < 6; ++i)
    {
        // TODO: DO this maybe
    }
}

void RockSaltOctahedraToggler::increment_leashes(const std::array<Coordinate, 6>& neighboring_vertex_ion_coords)
{
    for (const auto& vertex_ion_coord : neighboring_vertex_ion_coords)
    {
        auto vertex_ion_ix = this->coordinate_to_index(vertex_ion_coord);
        ++leashed_vertex_ions[vertex_ion_ix];
        assert(leashed_vertex_ions[vertex_ion_ix] <= 6);
    }
    return;
}

void RockSaltOctahedraToggler::reduce_leashes(const std::array<Coordinate, 6>& neighboring_vertex_ion_coords)
{
    for (const auto& vertex_ion_coord : neighboring_vertex_ion_coords)
    {
        auto vertex_ion_ix = this->coordinate_to_index(vertex_ion_coord);
        --leashed_vertex_ions[vertex_ion_ix];
        assert(leashed_vertex_ions[vertex_ion_ix] >= 0);
    }
    return;
}

void RockSaltOctahedraToggler::print(std::ostream& out_stream) const
{
    this->commit();
    Simplicity::print_poscar(this->rocksalt_struc, out_stream);
    return;
}

void RockSaltOctahedraToggler::commit() const
{
    this->commit_central_ions();
    this->commit_vertex_ions();

    return;
}

void RockSaltOctahedraToggler::commit_central_ions() const
{
    for (const auto& ix_is_on : this->central_ion_is_on)
    {
        auto ix = ix_is_on.first;
        auto is_on = ix_is_on.second;

        if (is_on)
        {
            ::set_site_occupant(&this->rocksalt_struc.basis[ix], this->central_ion_name);
        }

        else
        {
            ::set_site_occupant(&this->rocksalt_struc.basis[ix], "Va");
        }
    }
    return;
}

void RockSaltOctahedraToggler::commit_vertex_ions() const
{
    for (const auto& ix_count : this->leashed_vertex_ions)
    {
        auto ix = ix_count.first;
        auto count = ix_count.second;

        if (count > 0)
        {
            ::set_site_occupant(&this->rocksalt_struc.basis[ix], this->vertex_ion_name);
        }

        else
        {
            ::set_site_occupant(&this->rocksalt_struc.basis[ix], "Va");
        }
    }
    return;
}

RockSaltOctahedraToggler::index RockSaltOctahedraToggler::coordinate_to_index(const Coordinate& coordinate) const
{
    // Go through the basis of the structure
    // and find out which basis index the
    // given coordinate corresponds to
    const auto& basis = this->rocksalt_struc.basis;
    for (int ix = 0; ix < basis.size(); ++ix)
    {
        if (static_cast<CASM::Coordinate>(basis[ix]) == coordinate)
        {
            return ix;
        }
    }

    throw UtilExcept::IncompatibleCoordinate();
}

RockSaltOctahedraToggler::Coordinate RockSaltOctahedraToggler::index_to_coordinate(index coordinate_index) const
{
    // Go through the basis of the structure
    // and find out which coordinate the
    // given index corresponds to
    CASM::Site coord = rocksalt_struc.basis[coordinate_index];
    return coord;
    // TODO: explicitly cast for clairty?
}

RockSaltOctahedraToggler::Structure
RockSaltOctahedraToggler::primitive_structure(std::pair<std::string, std::string> species_names)
{
    const auto& central_ion_name = species_names.first;
    const auto& vertex_ion_name = species_names.second;

    auto nn_distance = RockSaltOctahedraToggler::nearest_neighbor_distance();
    auto lat = CASM::Lattice::fcc();
    auto scaled_lat = lat.scaled_lattice(nn_distance);

    CASM::Array<CASM::Site> basis;

    CASM::Coordinate pos_central(0, 0, 0, scaled_lat, CASM::FRAC);
    CASM::Coordinate pos_vertex(0.5, 0.5, 0.5, scaled_lat, CASM::FRAC);

    auto central_site = Extend::multi_atomic_site(pos_central, std::vector<std::string>{central_ion_name, "Va"});
    auto vertex_site = Extend::multi_atomic_site(pos_vertex, std::vector<std::string>{vertex_ion_name, "Va"});

    basis.push_back(central_site);
    basis.push_back(vertex_site);

    CASM::Structure primitive_structure(scaled_lat);
    primitive_structure.set_basis(basis);

    return primitive_structure;
}

RockSaltOctahedraToggler::Structure conventional_structure(std::pair<std::string, std::string> species_names,
                                                           std::string central_specie)
{
    // get the primitive,
    // apply trasnformation matrix
}

RockSaltOctahedraToggler RockSaltOctahedraToggler::relative_to_primitive(
    const Eigen::Matrix3i trans_mat, std::pair<std::string, std::string> species_names, bool central_is_second)
{
    if (central_is_second)
    {
        std::swap(species_names.first, species_names.second);
    }

    auto prim_struc = primitive_structure(species_names);
    auto super_struc = Simplicity::make_super_structure(prim_struc, trans_mat);
    auto super_lattice = super_struc.lattice();

    return RockSaltOctahedraToggler(std::move(super_struc), species_names.first, species_names.second,
                                    initialized_nearest_neighbor_deltas(super_lattice),
                                    initialized_central_ion_is_on(super_struc, species_names.first),
                                    initialized_leashed_vertex_ions(super_struc, species_names.second));
}

RockSaltOctahedraToggler RockSaltOctahedraToggler::relative_to_conventional(
    const Eigen::Matrix3i trans_mat, std::pair<std::string, std::string> species_names, bool central_is_second)
{
    //
    //
    //
}

RockSaltOctahedraToggler::RockSaltOctahedraToggler(Structure&& init_struc, std::string init_central_name,
                                                   std::string init_vertex_name,
                                                   std::array<Coordinate, 6>&& init_nn_deltas,
                                                   std::unordered_map<index, bool>&& init_central_is_on,
                                                   std::unordered_map<index, int>&& init_leashes)
    : rocksalt_struc(init_struc),
      central_ion_name(init_central_name),
      vertex_ion_name(init_vertex_name),
      nearest_neighbor_deltas(init_nn_deltas),
      central_ion_is_on(init_central_is_on),
      leashed_vertex_ions(init_leashes)
{
}

double RockSaltOctahedraToggler::nearest_neighbor_distance() { return 0.5; }

std::array<RockSaltOctahedraToggler::Coordinate, 6>
RockSaltOctahedraToggler::initialized_nearest_neighbor_deltas(const Lattice& rocksalt_lattice)
{
    // set of coordinates that take you from the central coordinate (center of octahedron)
    // to the six nearest nearest neighbor sites
    double d = RockSaltOctahedraToggler::nearest_neighbor_distance();
    // need a lattice for the coordinate type
    // need to CART type for coordinate, but this syntax is giving compiling errors...
    auto mode = CASM::CART;
    std::array<RockSaltOctahedraToggler::Coordinate, 6> deltas = {
        CASM::Coordinate(d, 0.0, 0.0, rocksalt_lattice, mode),
        CASM::Coordinate(0.0, d, 0.0, rocksalt_lattice, mode),
        CASM::Coordinate(0.0, 0.0, d, rocksalt_lattice, mode),
        CASM::Coordinate(-1 * d, 0.0, 0.0, rocksalt_lattice, mode),
        CASM::Coordinate(0.0, -1 * d, 0.0, rocksalt_lattice, mode),
        CASM::Coordinate(0.0, 0.0, -1 * d, rocksalt_lattice, mode)};

    return deltas;
}

std::unordered_map<RockSaltOctahedraToggler::index, bool>
RockSaltOctahedraToggler::initialized_central_ion_is_on(const Structure& init_struc, std::string central_name)
{
    std::unordered_map<RockSaltOctahedraToggler::index, bool> initialized_map;
    // TODO: range based loop
    for (int i = 0; i < init_struc.basis.size(); i++)
    {
        if (init_struc.basis[i].contains(central_name))
        {
            initialized_map[i] = false;
        }
    }

    return initialized_map;
}

std::unordered_map<RockSaltOctahedraToggler::index, int>
RockSaltOctahedraToggler::initialized_leashed_vertex_ions(const Structure& init_struc, std::string vertex_name)
{
    // I don't like this
    std::unordered_map<RockSaltOctahedraToggler::index, int> initialized_map;
    for (int i = 0; i < init_struc.basis.size(); i++)
    {
        if (init_struc.basis[i].contains(vertex_name)) // if basis is one of the vertex atoms
        {
            initialized_map[i] = 0;
        }
    }
    return initialized_map;
}

} // namespace SpecializedEnumeration
