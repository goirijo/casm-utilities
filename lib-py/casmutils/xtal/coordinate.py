from . import _xtal

Equals=_xtal.CoordinateEquals_f

class _Coordinate(_xtal.Coordinate):

    """Base class for both mutable and immutable Coordinate classes.
    Defines the functions that should be common for both."""

    def cart(self):
        """Return the Cartesian values of the coordinate
        Returns
        -------
        np.array

        """
        return self._cart_const()

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
        return self._frac_const(lat)

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
        self._equals=method(self, *args)

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
        return cls(_xtal.Coordinate.from_fractional(coords, lat)._cart_const())


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
        return self._equals(other)

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

class Coordinate(_Coordinate):

    """Immutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

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
        return self._bring_within_const(lat)

    def __add__(self, other):
        """Adds the "other" value to the Coordinate instance

        Parameters
        ----------
        other : Coordinate

        Returns
        -------
        Coordinate

        """
        return Coordinate(_xtal.Coordinate.__add__(self,other)._cart_const())

    def __iadd__(self, other):
        """Overloading the += operator defined the parent
        class to do nothing

        Parameters
        ----------
        other : Coordinate

        Returns
        -------
        TODO

        """
        pass


class MutableCoordinate(_Coordinate):

    """Mutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

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
        self._bring_within(lat)
        return

    def __add__(self, other):
        """Adds the "other" value to the
        MutableCoordinate instance

        Parameters
        ----------
        other : MutableCoordinate

        Returns
        -------
        MutableCoordinate

        """
        return MutableCoordinate(_xtal.Coordinate.__add__(self,other)._cart_const())

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
        return MutableCoordinate(_xtal.Coordinate.__iadd__(self,other)._cart_const())
