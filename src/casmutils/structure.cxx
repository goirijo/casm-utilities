#include "casmutils/structure.hpp"
#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <fstream>

namespace Rewrap
{
Structure::Structure(const CASM::Structure& init_struc) : CASM::Structure(init_struc) {}
Structure::Structure(Rewrap::fs::path& filename) : CASM::Structure(filename) {}

Structure Structure::from_poscar(const fs::path& poscar_path)
{
    return Rewrap::Structure(CASM::Structure(poscar_path));
}

bool Structure::is_primitive() const { return CASM::Structure::is_primitive(); }

Structure Structure::primitive() const { return Simplicity::make_primitive(*this); }
} // namespace Rewrap

namespace Simplicity
{
Rewrap::Structure make_primitive(const Rewrap::Structure& input)
{
    const CASM::Structure& casted_input(input);
    CASM::Structure true_prim;
    bool is_prim = casted_input.is_primitive(true_prim);
    return true_prim;
}

Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli)
{
    CASM::Structure niggli = non_niggli;
    CASM::Lattice lat_niggli = CASM::niggli(non_niggli.lattice(), CASM::TOL);
    niggli.set_lattice(lat_niggli, CASM::CART);
    return niggli;
}

void make_niggli(Rewrap::Structure* non_niggli)
{
    CASM::Lattice lat_niggli = CASM::niggli(non_niggli->lattice(), CASM::TOL);
    non_niggli->set_lattice(lat_niggli, CASM::CART);
    return;
}

void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(printable);
    p.sort();
    p.print(outstream);
    return;
}

void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}
} // namespace Simplicity
