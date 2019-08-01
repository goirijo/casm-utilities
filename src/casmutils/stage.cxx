#include "casmutils/stage.hpp"
#include "casmutils/exceptions.hpp"

namespace SpecializedEnumeration
{

void RockSaltOctahedraToggler::activate(const Coordinate& central_coord)
{
    auto central_coord_index = this->coordinate_to_index(central_coord);

    if (!this->is_cation_index(central_coord_index))
    {
        throw UtilExcept::IncompatibleCoordinate();
    }

    if (cation_is_on.at(central_coord_index))
    {
        return;
    }

    this->turn_cation_on(central_coord_index);

    auto nearest_anions = this->nearest_site_coordinates(central_coord);
    this->increment_leashes(nearest_anions);

    return;
}


void RockSaltOctahedraToggler::deactivate(const Coordinate& central_coord)
{
	auto central_coord_index = this->coordinate_to_index(central_coord);
        if (!this->is_cation_index(central_coord_index))
        {
             throw UtilExcept::IncompatibleCoordinate();
        }

        if (!cation_is_on.at(central_coord_index))
        {
             return;
        }
       
        this->turn_cation_off(central_coord_index);
        auto nearest_anions = this->nearest_site_coordinates(central_coord);
        this->reduce_leashes(nearest_anions);

       return;
}

bool RockSaltOctahedraToggler::is_cation_index(index ix) const { return cation_is_on.find(ix) != cation_is_on.end(); }

void RockSaltOctahedraToggler::turn_cation_on(index ix) { cation_is_on[ix] = true; }

void RockSaltOctahedraToggler::turn_cation_off(index ix) { cation_is_on[ix] = false; }

std::array<RockSaltOctahedraToggler::Coordinate, 6>
RockSaltOctahedraToggler::nearest_site_coordinates(const Coordinate& central_coord) const
{
    std::array<Coordinate, 6> nearest_coords(central_coord);
    for (int i = 0; i < 6; ++i)
    {
        nearest_coords[i] += nearest_neighbor_deltas[i];
    }

    return nearest_coords;
}

void RockSaltOctahedraToggler::increment_leashes(const std::array<Coordinate, 6>& neighboring_anion_coords)
{
    for(const auto& anion_coord : neighboring_anion_coords)
    {
        auto anion_ix=this->coordinate_to_index(anion_coord);
        ++leashed_anions[anion_ix];
    }
    return;
}

void RockSaltOctahedraToggler::reduce_leashes(const std::array<Coordinate, 6>& neighboring_anion_coords)
{
    for(const auto& anion_coord : neighboring_anion_coords)
    {
        auto anion_ix=this->coordinate_to_index(anion_coord);
        --leashed_anions[anion_ix];
    }
    return;
}


} // namespace SpecializedEnumeration
