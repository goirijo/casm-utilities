from . import _xtal
from . import globaldef

class Equals:

    """A Lattice compare method which returns true
    if all the lattice vectors of "reference" lattice
    are within a specified tolerance of the lattice vectors
    of "other" lattice"""

    def __init__(self, ref_lattice, tol):
        """Construct Equals from Lattice or MutableLattice &
        a given tolerance

        Parameters
        ----------
        ref_site : Lattice
        tol : double

        """
        self._LatticeEquals_f = _xtal.LatticeEquals_f(ref_lattice, tol)

    def __call__(self, other):
        """
        Parameters
        ----------
        other : Lattice

        Returns
        -------
        bool

        """
        return self._LatticeEquals_f(other)

class Lattice(_xtal.Lattice):

    """A Lattice class. Defined as the unit cells that
    define the unit cell"""

    def __init__(self, a, b, c):
        """
        Parameters
        ----------
        a : np.array
        b : np.array
        c : np.array

        """
        super().__init__(a,b,c)

    @classmethod
    def _from_pybind(cls, pybind_value):
        """Returns a constructed Lattice from
        a given _xtal.Lattice value

        Parameters
        ----------
        value : _xtal.Lattice

        Returns
        -------
        Lattice

        """
        value = cls(pybind_value.a(), pybind_value.b(), pybind_value.c())
        return value

    def a(self):
        """Returns 1st lattice vector

        Returns
        -------
        np.array

        """
        return super().a()

    def b(self):
        """Returns 2nd lattice vector

        Returns
        -------
        np.array

        """
        return super().b()

    def c(self):
        """Returns 3rd lattice vector

        Returns
        -------
        np.array

        """
        return super().c()

    def volume(self):
        """Returns volume of the lattice

        Returns
        -------
        double

        """
        return super().volume()

    def column_vector_matrix(self):
        """Returns the lattice in column vector
        format

        Returns
        -------
        np.array

        """
        return super().column_vector_matrix()

    def row_vector_matrix(self):
        """Returns the lattice in row vector
        format

        Returns
        -------
        np.array

        """
        return super().row_vector_matrix()

    def set_compare_method(self, method, *args):
        """Determines what strategy to use for comparing
        Lattices

        Parameters
        ----------
        method : Functor class that performs the computation
        *args : Arguments needed to construct the Functor

        """
        self._equals = method(self, *args)

    def __eq__(self, other):
        """Passes the "other" to the cuurent compare functor
        and compares it to the self and returns the evaluation

        Parameters
        ----------
        other : _Lattice

        Returns
        -------
        bool

        """
        if hasattr(self, '_equals') is False:
            self.set_compare_method(Equals, globaldef.tol)

        return self._equals(other)

    def __ne__(self, other):
        """Passes the "other" to the current compare functor
        and compares to the self and returns the opposite of the
        evaluation

        Parameters
        ----------
        other : _Lattice

        Returns
        -------
        bool

        """
        return not self == other

    def __str__(self):
        """Returns column vector matrix as a printable
        string

        Returns
        -------
        string

        """
        return self.__str__()

