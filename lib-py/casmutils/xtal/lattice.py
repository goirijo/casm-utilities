from . import _xtal
from . import globaldef

class Equals:

    """A wrapper class for _xtal.LatticeEquals_f"""

    def __init__(self, ref_lattice, tol):
        """Construct Equals from Lattice or MutableLattice &
        a given tolerance

        Parameters
        ----------
        ref_site : Lattice or MutableLattice
        tol : double

        """
        self._LatticeEquals_f = _xtal.LatticeEquals_f(ref_lattice._pybind_value, tol)

    def __call__(self, other):
        """Overloading () operator

        Parameters
        ----------
        other : Lattice or MutableLattice

        Returns
        -------
        bool

        """
        return self._LatticeEquals_f(other._pybind_value)

class _Lattice():

    """Base class for both mutable and immutable Lattice classes.
    Defines the functions that should be common for both."""

    def __init__(self, a, b, c):
        """create an instance of _xtal.Lattice as a container
        to access it's member function

        Parameters
        ----------
        a : np.array
        b : np.array
        c : np.array

        """
        if a is _xtal.Lattice and b is None and c is None:
            self._pybind_value = None

        else:
            self._pybind_value = _xtal.Lattice(a,b,c)

    @classmethod
    def _from_pybind(cls, py_bind_value):
        """Returns a constructed _Lattice from
        a given _xtal.Lattice value

        Parameters
        ----------
        value : _xtal.Lattice

        Returns
        -------
        _Lattice

        """
        value = cls(_xtal.Lattice,None,None)
        value._pybind_value = py_bind_value
        return value

    def a(self):
        """Returns 1st lattice vector

        Returns
        -------
        np.array

        """
        return self._pybind_value._a_const()

    def b(self):
        """Returns 2nd lattice vector

        Returns
        -------
        np.array

        """
        return self._pybind_value._b_const()

    def c(self):
        """Returns 3rd lattice vector

        Returns
        -------
        np.array

        """
        return self._pybind_value._c_const()

    def volume(self):
        """Returns volume of the lattice

        Returns
        -------
        double

        """
        return self._pybind_value._volume_const()

    def column_vector_matrix(self):
        """Returns the lattice in column vector
        format

        Returns
        -------
        np.array

        """
        return self._pybind_value._col_vec_mat_const()

    def row_vector_matrix(self):
        """Returns the lattice in row vector
        format

        Returns
        -------
        np.array

        """
        return self._pybind_value._row_vec_mat_const()

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

    def __getitem__(self, index):
        """Returns the lattice vector at the given
        index accessed through []

        Paremeters
        ----------
        index : int

        Returns
        -------
        np.array

        """
        return self._pybind_value[index]

    def __str__(self):
        """Returns the string to print

        Returns
        -------
        string

        """
        return self._pybind_value.__str__()

class Lattice(_Lattice):

    """Immutable Lattice class. Defined as three
    vectors that define the unit cell."""

    def __init__(self, a, b, c):
        """Constructor inheriting from
        _Lattice

        Paremeters
        ----------
        a : np.array
        b : np.array
        c : np.array

        """
        super().__init__(a, b, c)

class MutableLattice(_Lattice):

    """Mutable Lattice class. Defined as three
    vectors that define the unit cell"""

    def __init__(self, a, b, c):
        """Constructor inheriting from
        _Lattice

        Parameters
        ----------
        a : np.array
        b : np.array
        c : np.array

        """
        super().__init__(a, b, c)

