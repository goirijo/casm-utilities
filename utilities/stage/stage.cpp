#include "casmutils/stage.hpp"
#include "casmutils/structure_tools.hpp"
#include <iostream>

int main()
{
    Eigen::Matrix3i trans_mat;
    trans_mat<<-2,2,2,2,-2,2,2,2,-2;
    std::cout<<trans_mat<<std::endl;

    std::pair<std::string,std::string> species=std::make_pair("Na","Cl");

    auto struc=SpecializedEnumeration::RockSaltOctahedraToggler::primitive_structure(species);
    Simplicity::print_poscar(struc,std::cout);

    auto rocksalt=SpecializedEnumeration::RockSaltOctahedraToggler::relative_to_primitive(trans_mat, species, false);

    rocksalt.activate_all();
    rocksalt.deactivate(Rewrap::Coordinate(0,0,0));

    rocksalt.print(std::cout);

    return 0;
}
