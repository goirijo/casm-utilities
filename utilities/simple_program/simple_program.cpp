#include <iostream>
#include <fstream>
#include <casm/CASM_global_definitions.hh>
#include "casm/crystallography/Structure.hh"
#include <casm/casm_io/VaspIO.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/external/Eigen/Core>
#include <casm/crystallography/SupercellEnumerator.hh>

int main()
{
    /**
     * This is where you get to shine. Place all your code in here
     * and use the infrastructure of the repository to compile
     * and install.
     *
     * If you're interested in working with something more complicated
     * than a single main.cpp file is convenient for, check out
     * casm-utilities!
     * https://github.com/goirijo/casm-utilities
     */

    std::cout<<"I am a less sad executable that can do something"<<std::endl;
    CASM::print_splash(std::cout);

    std::cout<<"L E T ' S   D O   T H I S"<<std::endl;

    std::cout<<"- 1 - Read and write POSCAR"<<std::endl;
    
    // initialize structure from POSCAR file
    CASM::Structure my_structure("/home/julija/programming/casm-derived-template/src/test/POSCAR_1");
    
    // print POSCAR from structure object
    CASM::VaspIO::PrintPOSCAR my_printer(my_structure);
    
    // print to terminal
    my_printer.print(std::cout);
    
    // print to file
    std::ofstream new_poscar_file("/home/julija/programming/casm-derived-template/src/test/POSCAR_2");
    my_printer.print(new_poscar_file);
    new_poscar_file.close();

    std::cout<<"- 2 - Find prim"<<std::endl;
    
    // make PrimClex object
    CASM::PrimClex my_primclex(my_structure);
    CASM::Structure my_prim = my_primclex.get_prim();


    //CASM::ScelEnumProps my_scelenumprops(1,1);
    //my_primclex.generate_supercells(my_scelenumprops); // tricky
    //my_primclex.print_supercells(std::cout);

    // find prim
    //CASM::Structure my_prim = my_primclex.prim();
    CASM::VaspIO::PrintPOSCAR my_printer_prim(my_prim);
    my_printer_prim.print(std::cout);
    
    std::cout << "- 3 - Find sym group" << std::endl;
    
    auto my_fact_group = my_structure.factor_group();
    //my_fact_group.print(std::cout, CASM::FRAC);

    auto my_pt_group = my_structure.point_group();
    //my_pt_group.print(std::cout, CASM::FRAC);
    
    std::cout << "- 4 - Make Super" << std::endl;
    
    auto my_lat = my_structure.lattice();
    CASM::fs::ifstream TM_file("/home/julija/programming/casm-derived-template/src/test/TM");
    Eigen::Matrix3d TM;
    TM_file >> TM;
    TM_file.close();
    //TM << 2, 0, 0, 0, 2, 0, 0, 0, 2;
    auto &new_lat = my_lat.lat_column_mat()*TM;
    //std::cout << new_lat.col(0) << std::endl;
    auto super_lattice = CASM::Lattice(new_lat);
    //auto super_lattice = CASM::Lattice(new_lat.col(0), new_lat.col(1), new_lat.col(2));
    CASM::Structure super_structure = my_structure.create_superstruc(super_lattice);

    CASM::VaspIO::PrintPOSCAR my_printer_super(super_structure);
    my_printer_super.print(std::cout);

    return 0;
}
