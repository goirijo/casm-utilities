from . import _xtal
from .lattice import Lattice, MutableLattice
from .site import Site, MutableSite

class _Structure():

    """Base class for both mutable and immutable
    Structure classes. Defines function that are
    common to both"""

    def __init__(self, lattice, basis):
        """creates an instance of _xtal.Structure
        as a container to access it's member functions

        Parameters
        ----------
        lattice : Lattice or MutableLattice
        basis : list

        """
        if lattice is _xtal.Structure and basis is None:
            self._pybind_value = None

        else:
            py_bind_basis = []
            for site in basis:
                py_bind_basis.append(site._pybind_value)

            self._pybind_value = _xtal.Structure(lattice._pybind_value, py_bind_basis)

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
        Lattice or MutableLattice

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
        list

        """
        py_bind_basis = self._pybind_value._basis_sites_const()
        basis = []

        for py_bind_site in py_bind_basis:
            basis.append(Site._from_pybind(py_bind_site))

        return basis

    def __str__(self):
        """Prints the current structure instance

        Returns
        -------
        string

        """
        return self._pybind_value.__str__()

class Structure(_Structure):

    """Immutable structure defined by a Lattice and a list of basis sites"""

    def __init__(self, lattice, basis):
        """Inheriting the constructor from _Structure

        Paremeters
        ----------
        lattice : Lattice
        basis : list

        """
        super().__init__(lattice, basis)

    def set_lattice(self, new_lattice, coord_type):
        """Changes the lattice to the provided new lattice
        and updates the basis sites and returns a copy of
        the new structure

        Parameters
        ----------
        new_lattice : Lattice
        coord_type : string

        Returns
        -------
        Structure

        """
        return self._from_pybind(self._pybind_value._set_lattice_const(new_lattice._pybind_value, coord_type))

class MutableStructure(_Structure):

    """Mutable structure defined by a Lattice and a list of basis sites"""

    def __init__(self, lattice, basis):
        """Inheriting the constructor from _Structure

        Paremeters
        ----------
        lattice : Lattice or MutableLattice
        basis : list

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
        and updates the basis sites

        Parameters
        ----------
        new_lattice : Lattice or MutableLattice
        coord_type : string

        Returns
        -------
        TODO

        """
        self._pybind_value._set_lattice(new_lattice._pybind_value,coord_type)
        return

