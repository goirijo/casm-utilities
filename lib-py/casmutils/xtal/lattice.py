from . import _xtal

class _Lattice(_xtal.Lattice):

    """Base class for both mutable and immutable Lattice classes.
    Defines the functions that should be common for both."""

    pass

class Lattice(_Lattice):

    """Immutable Lattice class. Defined as three
    vectors that define the unit cell."""

    pass
