from . import _xtal

class Coordinate(_xtal.Coordinate):

    """Immutable Coordinate class. Defined as the Cartesian
    coodrinates, can handle opperations related to lattice
    periodicity."""

    def cart(self):
        """Return the Cartesian values of the coordinate
        Returns
        -------
        np.array

        """
        return self._cart_const()
