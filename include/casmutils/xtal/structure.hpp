#ifndef UTILS_STRUCTURE_HH
#define UTILS_STRUCTURE_HH

#include <casm/crystallography/BasicStructure.hh>
#include <casm/crystallography/SimpleStructure.hh>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <iostream>

namespace rewrap
{

/**
 * Describes a current, fixed state of a crystal. Composed of a lattice,
 * and collection of basis atoms (Sites). Each site contains the position of
 * the species, and a label for the species.
 */

class Structure
{
public:
    Structure() = delete;

    /// Construct by providing a path to a POSCAR like file
    static Structure from_poscar(const fs::path& poscar_path);

    /// Construct from parent CASM class
    Structure(const CASM::xtal::SimpleStructure& init_struc);
    Structure(const CASM::xtal::BasicStructure& init_struc);

    /// Construct with a lattice and list of sites (basis)
    Structure(const rewrap::Lattice& init_lat, const std::vector<Site>& init_basis);

    /// Returns a copy of the current lattice of the structure
    const Lattice& lattice() const;

    /// Give the structure a new lattice, and either keep the Cartesian, or fractional coordinates of the basis
    void set_lattice(const Lattice& new_lattice, COORD_TYPE mode);

    /// Give the structure a new lattice, and either keep the Cartesian, or fractional coordinates of the basis
    [[nodiscard]] Structure set_lattice(const Lattice& new_lattice, COORD_TYPE mode) const;

    /// Return a copy of all the basis sites
    const std::vector<Site>& basis_sites() const;

    /// Retreive the CASM implementations of *this
    template <typename CASMType> const CASMType& __get() const;

private:
    ///  updates the appropriate members using
    /// the template type given as a reference
    template <typename StructureRepresentation> void _update_using();
    ///  updates the casm_simplestructure member
    /// using casm_basicstructure as a reference
    void _update_simple_from_basic();
    ///  updates the casm_basicstructure member
    /// using casm_simplestructure as a reference
    void _update_basic_from_simple();

    ///  updates the lattice and basis member
    /// using casm_simplestructure as a reference
    void _update_internals_from_simple();

    ///  updates the lattice and basis member
    /// using casm_basicstructure as a reference
    void _update_internals_from_basic();

    ///  updates casm_basicstructure
    /// using the lattice and basis member as a reference
    void _update_basic_from_internals();

    ///  updates casm_simplestructure
    /// using the lattice and basis member as a reference
    void _update_simple_from_internals();

    /// CASM::SimpleStructure representation
    CASM::xtal::SimpleStructure casm_simplestructure;

    /// CASM::SimpleStructure representation
    CASM::xtal::BasicStructure casm_basicstructure;

    /// rewrap representation of the basis
    std::vector<Site> basis;

    /// rewrap representation of the lattice
    Lattice structure_lattice;
};
} // namespace rewrap

namespace casmutils
{
namespace xtal
{
using rewrap::CART;
using rewrap::COORD_TYPE;
using rewrap::FRAC;
using rewrap::Structure;
} // namespace xtal
} // namespace casmutils

#endif
