#include "casmutils/stage.hpp"
#include "casmutils/structure_tools.hpp"
#include <iostream>

int main()
{
    Eigen::Matrix3i trans_mat;
    trans_mat<<-1,3,4,1,-3,3,1,4,-2;
    std::cout<<trans_mat<<std::endl;

    std::pair<std::string,std::string> species=std::make_pair("Na","Cl");

    auto struc=SpecializedEnumeration::RockSaltOctahedraToggler::primitive_structure(species);
    Simplicity::print_poscar(struc,std::cout);

    auto rocksalt=SpecializedEnumeration::RockSaltOctahedraToggler::relative_to_primitive(trans_mat, species, false);

    auto central_coords=rocksalt.all_octahedron_center_coordinates();
    for(const auto ix_coord : central_coords)
    {
        std::cout<<ix_coord.first<<std::endl;
    }

    rocksalt.activate_all();

    /* for(int i=0; i<10; ++i) */
    /* { */
    /*     rocksalt.deactivate(i); */
    /* } */

    /* rocksalt.toggle_all(); */
    rocksalt.print(std::cout);

    return 0;
}
