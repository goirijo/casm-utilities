from . import _xtal

Equals=_xtal.CoordinateEquals_f

class _Coordinate:

    """Base class for both mutable and immutable Coordinate classes.
    Defines the functions that should be common for both."""

    def __init__(self, coord):
        """create an instance of _xtal.Coordinate
        as a container (_Coordinate object) which
        can be used to access the member functions
        of _xtal.Coordinate

        parameters
        ----------
        coord : np.array

        returns
        -------
        TODO

        """
        self._xtalCoordinate = _xtal.Coordinate(coord)

    def cart(self):
        """Return the Cartesian values of the coordinate

        Returns
        -------
        np.array

        """
        return self._xtalCoordinate._cart_const()

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
        return self._xtalCoordinate._frac_const(lat)

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
        Coordinate

        """
        return cls(_xtal.Coordinate.from_fractional(coords,lat)._cart_const())

    def set_compare_method(self, method, *args):
        """Determines what strategy should be used for comparison methods
        of Coordinates (e.g. compare Cartesian values within tolerance, or
        compare after bringing Coordinate within a unit cell).

        Parameters
        ----------
        method : Functor class that performs the evaluation
        *args : Arguments needed to construct method

        Returns
        -------
        TODO

        """
        self._equals=method(self._xtalCoordinate, *args)

    def __eq__(self, other):
        """Passes the "other" value to the current comparator
        stored in the Coordinate instance and returns the evaluation

        Parameters
        ----------
        other : Coordinate

        Returns
        -------
        bool

        """
        return self._equals(_xtal.Coordinate(other.cart()))

    def __ne__(self, other):
        """Passes the "other" value to the current comparator
        stored in the Coordinate instance and returns the
        opposite of the evaluation

        Parameters
        ----------
        other : Coordinate

        Returns
        -------
        bool

        """
        return not self==other

    def __add__(self, other):
        """Adds the "other" value to the Coordinate instance

        Parameters
        ----------
        other : Coordinate

        Returns
        -------
        Coordinate

        """
        return self.__class__(self._xtalCoordinate.__add__(_xtal.Coordinate(other.cart()))._cart_const())

class Coordinate(_Coordinate):

    """Immutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

    def __init__(self, coord):
        """Constructor inheriting from the
        parent _Coordinate

        parameters
        ----------
        coord : np.array

        returns
        -------
        TODO

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
        return self.__class__(self._xtalCoordinate._bring_within_const(lat)._cart_const())

class MutableCoordinate(_Coordinate):

    """Mutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

    def __init__(self, coord):
        """constructor that inherits from the
        parent _Coordinate

        parameters
        ----------
        coord : np.array

        returns
        -------
        TODO

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
        self._xtalCoordinate._bring_within(lat)
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
        return self.__class__(self._xtalCoordinate.__iadd__(_xtal.Coordinate(other.cart()))._cart_const())

