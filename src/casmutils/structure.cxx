#include <casmutils/structure.hpp>
#include <boost/filesystem.hpp>


namespace Rewrap
{
Structure::Structure(const CASM::Structure& init_struc) : CASM::Structure(init_struc) {}
Structure::Structure(Rewrap::fs::path& filename) : CASM::Structure(filename) {}

Structure Structure::from_poscar(const fs::path& poscar_path)
{
    return Rewrap::Structure(CASM::Structure(poscar_path));
}

bool Structure::is_primitive() const { return CASM::Structure::is_primitive(); }
} // namespace Rewrap

