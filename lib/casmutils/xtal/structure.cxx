#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/SimpleStructure.hh>
#include <casm/crystallography/SimpleStructureTools.hh>
#include <casm/global/definitions.hh>
#include <casmutils/definitions.hpp>
#include <casmutils/exceptions.hpp>
#include <casmutils/misc.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure.hpp>
#include <fstream>
namespace extend
{
} // namespace extend

namespace rewrap
{
Structure::Structure(const CASM::xtal::BasicStructure& init_struc)
    : structure_lattice(init_struc.lattice()), casm_basicstructure(init_struc)
{
    _update_using("basic");
}
Structure::Structure(const CASM::xtal::SimpleStructure& init_struc)
    : structure_lattice(init_struc.lat_column_mat), casm_simplestructure(init_struc)
{
    _update_using("simple");
}
Structure::Structure(const rewrap::Lattice& init_lat, const std::vector<rewrap::Site>& init_basis)
    : structure_lattice(init_lat), basis(init_basis)
{
    _update_using("internal");
}

Structure Structure::from_poscar(const fs::path& poscar_path)
{
    CASM::xtal::BasicStructure pos;
    if (!fs::exists(poscar_path))
    {
        throw except::UserInputMangle("Path does not exist");
    }
    std::ifstream infile(poscar_path);
    pos.read(infile);
    return casmutils::xtal::Structure(pos);
}

const Lattice& Structure::lattice() const { return this->structure_lattice; }

void Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode)
{
    this->casm_basicstructure.set_lattice(CASM::xtal::Lattice(new_lattice.column_vector_matrix()), mode);
    _update_using("basic");
    return;
}

[[nodiscard]] Structure Structure::set_lattice(const Lattice& new_lattice, COORD_TYPE mode) const
{
    casmutils::xtal::Structure new_struc = *this;
    new_struc.set_lattice(new_lattice, mode);
    return new_struc;
}
const std::vector<Site>& Structure::basis_sites() const { return this->basis; }

void Structure::_update_using(std::string reference)
{
    if (reference == "basic")
    {
        _update_simple_from_basic();
        _update_internals_from_basic();
    }
    else if (reference == "simple")
    {
        _update_basic_from_simple();
        _update_internals_from_simple();
    }
    else if (reference == "internal")
    {
        _update_basic_from_internals();
        _update_simple_from_internals();
    }
    else
    {
        throw except::UserInputMangle("Check the update reference you used. Not compatible.");
    }
}

void Structure::_update_simple_from_basic()
{
    this->casm_simplestructure = CASM::xtal::make_simple_structure(this->casm_basicstructure);
}

void Structure::_update_internals_from_basic()
{
    this->structure_lattice = casmutils::xtal::Lattice(this->casm_basicstructure.lattice());
    this->basis.clear();
    for (const auto& site : this->casm_basicstructure.basis())
    {
        this->basis.emplace_back(site, 0);
    }
}

void Structure::_update_basic_from_simple()
{
    /// Syncing Lattice
    this->casm_basicstructure =
        CASM::xtal::BasicStructure(CASM::xtal::Lattice(this->casm_simplestructure.lat_column_mat));
    /// Lattice has been synced
    auto& basic_basis = this->casm_basicstructure.set_basis();
    for (int index = 0; index < this->casm_simplestructure.n_atom(); index++)
    {
        const auto& info = this->casm_simplestructure.info(CASM::xtal::SimpleStructure::SpeciesMode::ATOM);
        Eigen::Vector3d raw_coord = info.coord(index);
        basic_basis.emplace_back(CASM::xtal::Coordinate(raw_coord, this->casm_basicstructure.lattice(), CASM::CART),
                                 info.names[index]);
    }
}

void Structure::_update_internals_from_simple()
{
    this->structure_lattice = casmutils::xtal::Lattice(casm_simplestructure.lat_column_mat);
    this->basis.clear();
    for (int index = 0; index < this->casm_simplestructure.n_atom(); index++)
    {
        const auto& info = this->casm_simplestructure.info(CASM::xtal::SimpleStructure::SpeciesMode::ATOM);
        Eigen::Vector3d raw_coord = info.coord(index);
        this->basis.emplace_back(casmutils::xtal::Coordinate(raw_coord), info.names[index]);
    }
}

void Structure::_update_basic_from_internals()
{
    /// Syncing Lattice
    this->casm_basicstructure =
        CASM::xtal::BasicStructure(CASM::xtal::Lattice(this->structure_lattice.column_vector_matrix()));
    /// Lattice has been synced
    auto& basic_basis = this->casm_basicstructure.set_basis();
    for (const auto& site : this->basis)
    {
        basic_basis.emplace_back(CASM::xtal::Coordinate(site.cart(), this->casm_basicstructure.lattice(), CASM::CART),
                                 site.label());
    }
}

void Structure::_update_simple_from_internals()
{
    /// This seems like a horrible way to do this but it's because SimpleStructure only has a
    // default constructor and a whole bunch of empty fields.
    // I also didn't want to call _update_basic_from_internals here.
    CASM::xtal::BasicStructure temp_basic(CASM::xtal::Lattice(this->structure_lattice.column_vector_matrix()));
    auto& basic_basis = temp_basic.set_basis();
    for (const auto& site : this->basis)
    {
        basic_basis.emplace_back(CASM::xtal::Coordinate(site.cart(), temp_basic.lattice(), CASM::CART), site.label());
    }
    this->casm_simplestructure = CASM::xtal::make_simple_structure(temp_basic);
}

/// Return *this as a CASM::BasicStructure
template <> const CASM::xtal::SimpleStructure& Structure::__get<CASM::xtal::SimpleStructure>() const
{
    return this->casm_simplestructure;
}

/// Return *this as a CASM::BasicStructure
template <> const CASM::xtal::BasicStructure& Structure::__get<CASM::xtal::BasicStructure>() const
{
    return this->casm_basicstructure;
}

} // namespace rewrap
