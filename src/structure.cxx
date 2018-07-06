#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <boost/filesystem.hpp>
#include <fstream>

namespace Simplicity
{
CASM::Structure make_primitive(const CASM::Structure& input)
{
    CASM::Structure true_prim;
    bool is_prim = input.is_primitive(true_prim);
    return true_prim;
}

/* CASM::Lattice canonical_equivalent_lattice(CASM::Structure& non_equiv) */
/* { */
/*     return CASM::canonical_equivalent_lattice(non_equiv.lattice(), non_equiv.point_group(), CASM::TOL); */
/* } */

CASM::Structure make_niggli(const CASM::Structure& non_niggli)
{
    CASM::Structure niggli = non_niggli;
    CASM::Lattice lat_niggli = CASM::niggli(non_niggli.lattice(), CASM::TOL);
    niggli.set_lattice(lat_niggli, CASM::CART);
    return niggli;
}

void print_poscar(const CASM::Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(printable);
    p.sort();
    p.print(outstream);
    return;
}

void print_poscar(const CASM::Structure& printable, const CASM::fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}
}
