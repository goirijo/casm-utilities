#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/rocksalttoggler.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>

namespace
{
// Please never use this outisde of the RockSalt context
void set_site_occupant(CASM::xtal::Site* mutating_site, std::string new_occ)
{
    throw except::NotImplemented();
    /* if (new_occ == "Va") */
    /* { */
    /*     mutating_site->set_occ_value(0); */
    /*     assert(mutating_site->allowed_occupants()[0] == "Va"); */
    /* } */

    /* else */
    /* { */
    /*     mutating_site->set_occ_value(1); */
    /* } */
}
} // namespace

namespace enumeration
{

void RockSaltOctahedraToggler::activate(index central_coord_index)
{
    if (!this->is_central_ion_index(central_coord_index))
    {
        throw except::IncompatibleCoordinate();
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

void RockSaltOctahedraToggler::activate_all()
{
    for (const auto& centralix_ison_pair : this->central_ion_is_on)
    {
        auto ix = centralix_ison_pair.first;
        auto is_on = centralix_ison_pair.second;
        if (!is_on)
        {
            this->activate(ix);
        }
    }
    return;
}

void RockSaltOctahedraToggler::activate(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);
    this->activate(central_coord_index);
    return;
}

void RockSaltOctahedraToggler::toggle(index central_coord_index)
{
    if (!this->is_central_ion_index(central_coord_index))
    {
        throw except::IncompatibleCoordinate();
    }

    if (central_ion_is_on.at(central_coord_index))
    {
        this->deactivate(central_coord_index);
    }

    else
    {
        this->activate(central_coord_index);
    }

    return;
}

void RockSaltOctahedraToggler::toggle(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);
    this->toggle(central_coord_index);
    return;
}

void RockSaltOctahedraToggler::toggle_all()
{
    for (const auto& centralix_ison_pair : this->central_ion_is_on)
    {
        auto ix = centralix_ison_pair.first;
        this->toggle(ix);
    }
    return;
}

void RockSaltOctahedraToggler::deactivate(index central_coordinate_index)
{
    if (!this->is_central_ion_index(central_coordinate_index))
    {
        throw except::IncompatibleCoordinate();
    }

    if (!central_ion_is_on.at(central_coordinate_index))
    {
        return;
    }

    this->turn_central_ion_off(central_coordinate_index);

    auto nearest_vertex_ions =
        this->nearest_neighbor_site_coordinates(this->index_to_coordinate(central_coordinate_index));
    this->reduce_leashes(nearest_vertex_ions);

    return;
}

void RockSaltOctahedraToggler::deactivate(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);
    this->deactivate(central_coord_index);
    return;
}

void RockSaltOctahedraToggler::deactivate_all()
{
    for (const auto& centralix_ison_pair : this->central_ion_is_on)
    {
        auto ix = centralix_ison_pair.first;
        auto is_on = centralix_ison_pair.second;
        if (is_on)
        {
            this->deactivate(ix);
        }
    }
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
    std::array<Coordinate, 6> nearest_coords{
        central_coord, central_coord, central_coord, central_coord, central_coord, central_coord};
    for (int i = 0; i < 6; ++i)
    {
        nearest_coords[i] += nearest_neighbor_deltas[i];
    }

    return nearest_coords;
}

std::array<RockSaltOctahedraToggler::index, 6>
RockSaltOctahedraToggler::nearest_neighbor_site_indexes(index central_coord_index) const
{
    throw except::NotImplemented();
    /* auto nearest_coordinates =
     * this->nearest_neighbor_site_coordinates(this->index_to_coordinate(central_coord_index)); */
    /* std::array<index, 6> nearest_coord_indexes; */

    /* for (int i = 0; i < 6; ++i) */
    /* { */
    /*     // TODO: DO this maybe */
    /* } */
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

    casmutils::xtal::print_poscar(this->rocksalt_struc, out_stream);
    return;
}

RockSaltOctahedraToggler::Structure RockSaltOctahedraToggler::structure() const
{
    this->commit();
    return this->rocksalt_struc;
}

void RockSaltOctahedraToggler::commit() const
{
    this->commit_central_ions();
    this->commit_vertex_ions();

    return;
}

void RockSaltOctahedraToggler::commit_central_ions() const
{
    throw except::NotImplemented();
    /* for (const auto& ix_is_on : this->central_ion_is_on) */
    /* { */
    /*     auto ix = ix_is_on.first; */
    /*     auto is_on = ix_is_on.second; */

    /*     if (is_on) */
    /*     { */
    /*         ::set_site_occupant(&this->rocksalt_struc.basis[ix], this->central_ion_name); */
    /*     } */

    /*     else */
    /*     { */
    /*         ::set_site_occupant(&this->rocksalt_struc.basis[ix], "Va"); */
    /*     } */
    /* } */
    return;
}

void RockSaltOctahedraToggler::commit_vertex_ions() const
{
    throw except::NotImplemented();
    /* for (const auto& ix_count : this->leashed_vertex_ions) */
    /* { */
    /*     auto ix = ix_count.first; */
    /*     auto count = ix_count.second; */

    /*     if (count > 0) */
    /*     { */
    /*         ::set_site_occupant(&this->rocksalt_struc.basis[ix], this->vertex_ion_name); */
    /*     } */

    /*     else */
    /*     { */
    /*         ::set_site_occupant(&this->rocksalt_struc.basis[ix], "Va"); */
    /*     } */
    /* } */
    return;
}

RockSaltOctahedraToggler::index RockSaltOctahedraToggler::coordinate_to_index(Coordinate coordinate) const
{
    throw except::NotImplemented();
    // TODO: Bring the coordinate within relative to the rocksalt structure lattice before making any comparisons
    //    coordinate.bring_within(this->rocksalt_struc.lattice());
    //    // Go through the basis of the structure
    //    // and find out which basis index the
    //    // given coordinate corresponds to
    //    auto basis = this->rocksalt_struc.basis_sites();
    //    for (int ix = 0; ix < basis.size(); ++ix)
    //    {
    //        if (casmutils::is_equal<casmutils::xtal::CoordinateEquals_f>(
    //                static_cast<Coordinate>(basis[ix]), coordinate, 1e-5))
    //        {
    //            return ix;
    //        }
    //    }
    //
    //    throw except::IncompatibleCoordinate();
}

RockSaltOctahedraToggler::Coordinate RockSaltOctahedraToggler::index_to_coordinate(index coordinate_index) const
{
    throw except::NotImplemented();
    /* // Go through the basis of the structure */
    /* // and find out which coordinate the */
    /* // given index corresponds to */
    /* const auto& coord = rocksalt_struc.basis[coordinate_index]; */
    /* return coord; */
    /* // TODO: explicitly cast for clairty? */
}

RockSaltOctahedraToggler::Structure
RockSaltOctahedraToggler::primitive_structure(std::pair<std::string, std::string> species_names,
                                              double init_nn_distance)
{
    throw except::NotImplemented();
    // const auto& central_ion_name = species_names.first;
    // const auto& vertex_ion_name = species_names.second;

    // Eigen::Matrix3d lat_mat;
    // lat_mat<<0,1,1,1,0,1,1,1,0;
    // lat_mat*=init_nn_distance;
    ///* auto scaled_lat = lat.scaled_lattice(2*nn_distance); */
    // Lattice scaled_lat(lat_mat);

    // std::vector<rewrap::Site> basis;

    // rewrap::Coordinate pos_central = rewrap::Coordinate::from_fractional(0, 0, 0, scaled_lat);
    // rewrap::Coordinate pos_vertex = rewrap::Coordinate::from_fractional(0.5, 0.5, 0.5, scaled_lat);

    // rewrap::Site central_site(pos_central, std::vector<std::string>{"Va", central_ion_name});
    // rewrap::Site vertex_site(pos_vertex, std::vector<std::string>{"Va", vertex_ion_name});

    // basis.push_back(central_site);
    // basis.push_back(vertex_site);

    // rewrap::Structure primitive_structure(scaled_lat, basis);
    // return primitive_structure;
}

RockSaltOctahedraToggler::Structure conventional_structure(std::pair<std::string, std::string> species_names,
                                                           std::string central_specie,
                                                           double init_nn_distance)
{
    throw except::NotImplemented();
    // get the primitive,
    // apply trasnformation matrix
}

RockSaltOctahedraToggler
RockSaltOctahedraToggler::relative_to_primitive(const Eigen::Matrix3i trans_mat,
                                                std::pair<std::string, std::string> species_names,
                                                bool central_is_second,
                                                double init_nn_distance)
{
    if (central_is_second)
    {
        std::swap(species_names.first, species_names.second);
    }

    auto prim_struc = primitive_structure(species_names, init_nn_distance);
    auto super_struc = casmutils::xtal::make_superstructure(prim_struc, trans_mat);

    return RockSaltOctahedraToggler(std::move(super_struc),
                                    species_names.first,
                                    species_names.second,
                                    initialized_nearest_neighbor_deltas(init_nn_distance),
                                    initialized_central_ion_is_on(super_struc, species_names.first),
                                    initialized_leashed_vertex_ions(super_struc, species_names.second),
                                    init_nn_distance);
}

RockSaltOctahedraToggler
RockSaltOctahedraToggler::relative_to_conventional(const Eigen::Matrix3i trans_mat,
                                                   std::pair<std::string, std::string> species_names,
                                                   bool central_is_second,
                                                   double init_nn_distance)
{
    throw except::NotImplemented();
}

RockSaltOctahedraToggler::RockSaltOctahedraToggler(Structure&& init_struc,
                                                   std::string init_central_name,
                                                   std::string init_vertex_name,
                                                   std::array<Coordinate, 6>&& init_nn_deltas,
                                                   std::unordered_map<index, bool>&& init_central_is_on,
                                                   std::unordered_map<index, int>&& init_leashes,
                                                   double init_nn_distance)
    : rocksalt_struc(init_struc),
      central_ion_name(init_central_name),
      vertex_ion_name(init_vertex_name),
      nearest_neighbor_deltas(init_nn_deltas),
      central_ion_is_on(init_central_is_on),
      leashed_vertex_ions(init_leashes),
      nn_distance(init_nn_distance)
{
}

double RockSaltOctahedraToggler::nearest_neighbor_distance() const { return this->nn_distance; }

std::array<RockSaltOctahedraToggler::Coordinate, 6>
RockSaltOctahedraToggler::initialized_nearest_neighbor_deltas(double init_nn_distance)
{
    auto d = init_nn_distance;
    // need a lattice for the coordinate type
    // need to CART type for coordinate, but this syntax is giving compiling errors...
    std::array<RockSaltOctahedraToggler::Coordinate, 6> deltas = {Coordinate(d, 0.0, 0.0),
                                                                  Coordinate(0.0, d, 0.0),
                                                                  Coordinate(0.0, 0.0, d),
                                                                  Coordinate(-1 * d, 0.0, 0.0),
                                                                  Coordinate(0.0, -1 * d, 0.0),
                                                                  Coordinate(0.0, 0.0, -1 * d)};

    return deltas;
}

std::unordered_map<RockSaltOctahedraToggler::index, bool>
RockSaltOctahedraToggler::initialized_central_ion_is_on(const Structure& init_struc, std::string central_name)
{
    throw except::NotImplemented();
    /* std::unordered_map<RockSaltOctahedraToggler::index, bool> initialized_map; */
    /* // TODO: range based loop */
    /* for (int i = 0; i < init_struc.basis.size(); i++) */
    /* { */
    /*     if (init_struc.basis[i].contains(central_name)) */
    /*     { */
    /*         initialized_map[i] = false; */
    /*     } */
    /* } */

    /* return initialized_map; */
}

std::unordered_map<RockSaltOctahedraToggler::index, int>
RockSaltOctahedraToggler::initialized_leashed_vertex_ions(const Structure& init_struc, std::string vertex_name)
{
    throw except::NotImplemented();
    /* // I don't like this */
    /* std::unordered_map<RockSaltOctahedraToggler::index, int> initialized_map; */
    /* for (int i = 0; i < init_struc.basis.size(); i++) */
    /* { */
    /*     if (init_struc.basis[i].contains(vertex_name)) // if basis is one of the vertex atoms */
    /*     { */
    /*         initialized_map[i] = 0; */
    /*     } */
    /* } */
    /* return initialized_map; */
}

std::vector<std::pair<RockSaltOctahedraToggler::index, RockSaltOctahedraToggler::Coordinate>>
RockSaltOctahedraToggler::all_octahedron_center_coordinates() const // Finish
{
    std::vector<std::pair<index, Coordinate>> index_coord_pair_list;
    for (const auto& octahedral_pairs : central_ion_is_on)
    {
        index ix = octahedral_pairs.first;
        auto coord = index_to_coordinate(ix);

        index_coord_pair_list.emplace_back(ix, coord);
    }

    return index_coord_pair_list;
}
} // namespace enumeration
