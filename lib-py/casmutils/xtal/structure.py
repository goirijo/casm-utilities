from . import _xtal
from .lattice import Lattice
from .site import Site, MutableSite

class _Structure():

    """Base class for both mutable and immutable
    Structure classes. Defines function that are
    common to both"""

    def __init__(self, lattice, basis):
        """
        Parameters
        ----------
        lattice : Lattice or MutableLattice
        basis : list[Sites]

        """
        if lattice is _xtal.Structure and basis is None:
            self._pybind_value = None

        else:
            py_bind_basis = [site._pybind_value for site in basis]
            self._pybind_value = _xtal.Structure(lattice, py_bind_basis)

    @classmethod
    def _from_pybind(cls, py_bind_value):
        """Returns a constructed _Structure from
        a given _xtal.Structure value

        Paremeters
        ----------
        py_bind_value : _xtal.Structure

        Returns
        -------
        _Structure

        """
        value = cls(_xtal.Structure, None)
        value._pybind_value = py_bind_value
        return value

    @classmethod
    def from_poscar(cls, poscar_path):
        """Reads POSCAR and returns
        constructed Structure

        Parameters
        ----------
        poscar_path : string

        Returns
        -------
        Structure or MutableStructure

        """
        py_bind_structure = _xtal.Structure._from_poscar(poscar_path)
        return cls._from_pybind(py_bind_structure)

    def lattice(self):
        """Returns the lattice of the
        current structure

        Returns
        -------
        Lattice

        """
        return Lattice._from_pybind(self._pybind_value._lattice_const())

    def to_poscar(self, filename):
        """Prints out the structure to a file
        in POSCAR format

        Parameters
        ----------
        filename : string

        """
        return self._pybind_value._to_poscar(filename)

    def basis_sites(self):
        """Returns the basis sites in the structure

        Returns
        -------
        list[Sites]

        """
        py_bind_basis = self._pybind_value._basis_sites_const()
        basis = [Site._from_pybind(py_bind_site) for py_bind_site in py_bind_basis]

        return basis

    def __str__(self):
        """Returns lattice and the list of basis sites as
        a printable string

        Returns
        -------
        string

        """
        return self._pybind_value.__str__()

class Structure(_Structure):

    """Immutable structure defined by a Lattice and a list of basis sites.
    Handles all the operations in a const way"""

    def __init__(self, lattice, basis):
        """
        Paremeters
        ----------
        lattice : Lattice
        basis : list[Sites]

        """
        super().__init__(lattice, basis)

    def set_lattice(self, new_lattice, coord_type):
        """Changes the lattice to the provided new lattice
        and updates the basis sites based on the coord_type provided
        (Fractional or Cartesian) and returns a copy of
        the new structure

        Parameters
        ----------
        new_lattice : Lattice
        coord_type : string

        Returns
        -------
        Structure

        """
        return self._from_pybind(self._pybind_value._set_lattice_const(new_lattice, coord_type))

class MutableStructure(_Structure):

    """Mutable structure defined by a Lattice and a list of basis sites.
    can handle non const operations. use only when you want to mutate
    the class itself."""

    def __init__(self, lattice, basis):
        """
        Paremeters
        ----------
        lattice : Lattice
        basis : list[Sites]

        """
        super().__init__(lattice, basis)

    def within(self):
        """Moves the basis sites within the lattice

        Returns
        -------
        None

        """
        self._pybind_value._within()
        return

    def set_lattice(self, new_lattice, coord_type):
        """Changes the lattice to the provided new lattice
        and updates the basis sites based on the coord_type
        (Fractional or Cartesian) provided

        Parameters
        ----------
        new_lattice : Lattice or MutableLattice
        coord_type : string

        Returns
        -------
        TODO

        """
        self._pybind_value._set_lattice(new_lattice,coord_type)
        return

