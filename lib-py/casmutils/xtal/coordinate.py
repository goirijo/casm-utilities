from . import _xtal
from . import globaldef
import numpy as np


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
        self.ref_coordinate = ref_coordinate
        self._CoordinateEquals_f = _xtal.CoordinateEquals_f(tol)

    def __call__(self, other):
        """
        Parameters
        ----------
        other : Coordinate or MutableCoordinate

        Returns
        -------
        bool

        """
        return self._CoordinateEquals_f(self.ref_coordinate.cart(),
                                        other.cart())


class _Coordinate:
    """Base class for both mutable and immutable Coordinate classes.
    Defines the functions that should be common for both."""
    def __init__(self, coord):
        """
        Parameters
        ----------
        coord : np.array

        """
        self.coord = coord

    def cart(self):
        """Returns the Cartesian values of the coordinate

        Returns
        -------
        np.array

        """
        return self.coord

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
        return _xtal.cartesian_to_fractional(self.coord, lat)

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
        cartesian_coords = _xtal.fractional_to_cartesian(coords, lat)
        return cls(cartesian_coords)

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
        return other.__class__(np.add(self.coord, other.cart()))

    def __str__(self):
        """Returns x,y,z values of the Coordinate as a printabe
        string

        Returns
        -------
        string

        """
        return self.coord.__str__()

    def __rmul__(self, CartOp):
        """Applying a given symop to Coordinate

        Parameters
        ----------
        CartOp : cu.sym.CartOp

        Returns
        -------
        Coordinate or MutableCoordinate

        """
        return self.__class__(CartOp * self.coord)


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
        return Coordinate(_xtal.bring_within_lattice(self.coord, lat))


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
        self.coord = _xtal.bring_within_lattice(self.coord, lat)
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
        self.coord = np.add(self.coord, other.cart())
        return self


def cartesian_to_fractional(cart_coords, lat):
    """Returs fractional coordinates of the given cartesian coordinates

    Parameters
    ----------
    cart_coords : np.array
    lat : cu.xtal.Lattice

    Returns
    -------
    np.array

    """
    return _xtal.cartesian_to_fractional(cart_coords, lat)


def fractional_to_cartesian(frac_coords, lat):
    """Returs fractional coordinates of the given cartesian coordinates

    Parameters
    ----------
    frac_coords: np.array
    lat : cu.xtal.Lattice

    Returns
    -------
    np.array

    """
    return _xtal.fractional_to_cartesian(frac_coords, lat)


def bring_within_lattice(cart_coords, lat):
    """Brings the given cartesian coordinates within the lattice

    Parameters
    ----------
    cart_coords : np.array
    lat : cu.xtal.Lattice

    Returns
    -------
    np.array

    """
    return _xtal.bring_within_lattice(cart_coords, lat)


def bring_within_wigner_seitz(cart_coords, lat):
    """Brings the given cartesian coordinates within the wigner seitz cell of the lattice

    Parameters
    ----------
    cart_coords : np.array
    lat : cu.xtal.Lattice

    Returns
    -------
    np.array

    """
    return _xtal.bring_within_wigner_seitz(cart_coords, lat)
