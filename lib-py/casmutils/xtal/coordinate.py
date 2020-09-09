from . import _xtal
from . import globaldef


class Equals:

    """A Coordinate compare method which returns true if
    each of the x,y,z values of the "reference" Coordinate
    are within a specified tolerance of x,y,z values
    of the "other" Coordinate or MutableCoordinate"""

    def __init__(self, ref_coordinate, tol):
        """Construct Equals from Coordinate
        or MutableCoordinate & a given tolerance

        Parameters
        ----------
        ref_coordinate : Coordinate or MutableCoordinate
        tol : double

        """
        self._CoordinateEquals_f = _xtal.CoordinateEquals_f(
            ref_coordinate._pybind_value, tol)

    def __call__(self, other):
        """
        Parameters
        ----------
        other : Coordinate or MutableCoordinate

        Returns
        -------
        bool

        """
        return self._CoordinateEquals_f(other._pybind_value)


class _Coordinate:

    """Base class for both mutable and immutable Coordinate classes.
    Defines the functions that should be common for both."""

    def __init__(self, coord):
        """
        Parameters
        ----------
        coord : np.array

        """
        if coord is _xtal.Coordinate:
            self._pybind_value = None

        else:
            self._pybind_value = _xtal.Coordinate(coord)

    def cart(self):
        """Returns the Cartesian values of the coordinate

        Returns
        -------
        np.array

        """
        return self._pybind_value._cart_const()

    def frac(self, lat):
        """Returns the fractional values of the coordinate
        relative to the given lattice

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        np.array

        """
        return self._pybind_value._frac_const(lat)

    @classmethod
    def _from_pybind(cls, py_bind_value):
        """Returns a constructed _Coordinate from
        a given _xtal.Coorindate value

        Paremeters
        ----------
        py_bind_value : _xtal.Coordinate

        Returns
        -------
        _Coordinate

        """
        value = cls(_xtal.Coordinate)
        value._pybind_value = py_bind_value
        return value

    @classmethod
    def from_fractional(cls, coords, lat):
        """Constructs a Coordinate from
        fractional coordinates

        Parameters
        ----------
        coords : np.array
        lat : Lattice

        Returns
        -------
        cls

        """
        py_binded = _xtal.Coordinate.from_fractional(coords, lat)
        return cls._from_pybind(py_binded)

    def set_compare_method(self, method, *args):
        """Determines what strategy should be used for comparison methods
        of Coordinates (e.g. compare Cartesian values within tolerance, or
        compare after bringing Coordinate within a unit cell).

        Parameters
        ----------
        method : Functor class that performs the evaluation
        *args : Arguments needed to construct method

        """
        self._equals = method(self, *args)

    def __eq__(self, other):
        """Passes the "other" value to the current comparator
        stored in the Coordinate instance and returns the evaluation

        Parameters
        ----------
        other : Coordinate or MutableCoordinate

        Returns
        -------
        bool

        """
        if hasattr(self, '_equals') is False:
            self.set_compare_method(Equals, globaldef.tol)

        return self._equals(other)

    def __ne__(self, other):
        """Passes the "other" value to the current comparator
        stored in the Coordinate instance and returns the
        opposite of the evaluation

        Parameters
        ----------
        other : Coordinate or MutablelCoordinate

        Returns
        -------
        bool

        """
        return not self == other

    def __add__(self, other):
        """Adds the "other" value to the Coordinate instance

        Parameters
        ----------
        other : Coordinate or MutableCoordinate

        Returns
        -------
        Coordinate or MutableCoordinate

        """
        py_binded = self._pybind_value + other._pybind_value
        return self._from_pybind(py_binded)

    def __str__(self):
        """Returns x,y,z values of the Coordinate as a printabe
        string

        Returns
        -------
        string

        """
        return self._pybind_value.__str__()

    def __rmul__(self, CartOp):
        """Applies the provided symmetry operation to the 
        Coordinate and returns the transformed coordinate

        Parameters
        ----------
        CartOp : cu.sym.CartOp

        Returns
        -------
        Coordinate

        """
        return self._from_pybind(CartOp * self._pybind_value)


class Coordinate(_Coordinate):

    """Immutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

    def __init__(self, coord):
        """
        Parameters
        ----------
        coord : np.array

        """
        super().__init__(coord)

    def bring_within(self, lat):
        """Returns the coordinate after applying lattice
        translations that bring it within the unit cell of the
        given lattice

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        Coordinate

        """
        py_binded = self._pybind_value._bring_within_const(lat)
        return self._from_pybind(py_binded)


class MutableCoordinate(_Coordinate):

    """Mutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

    def __init__(self, coord):
        """
        Parameters
        ----------
        coord : np.array

        """
        super().__init__(coord)

    def bring_within(self, lat):
        """Apply lattice translations to self
        that bring it within the unit cell of the
        given lattice

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        None

        """
        self._pybind_value._bring_within(lat)
        return

    def __iadd__(self, other):
        """Adds the "other" value to the current MutableCoordinate
        instance and makes it the cuurent instance

        Parameters
        ----------
        other : MutableCoordinate

        Returns
        -------
        MutableCoordinate

        """
        self._pybind_value += other._pybind_value
        return self
