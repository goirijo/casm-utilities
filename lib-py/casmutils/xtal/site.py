from . import _xtal
from . import globaldef


class Equals:

    """Definition of a site compare method which returns true if
    the coordinate and the atom type of the "reference" site matches
    with the coordinate and the atom type of the "other" site"""

    def __init__(self, ref_site, tol):
        """Construct Equals from Site
        or MutableSite and a given tolerance

        Parameters
        ----------
        ref_site : Site or MutableSite
        tol : double

        """
        self._SiteEquals_f = _xtal.SiteEquals_f(ref_site._pybind_value, tol)

    def __call__(self, other):
        """
        Parameters
        ----------
        other : Site or MutableSite

        Returns
        -------
        bool

        """
        return self._SiteEquals_f(other._pybind_value)


class _Site:

    """Base class for both mutable and immutable Site classes.
    Defines the functions that should be common to both"""

    def __init__(self, coord, label):
        """
        Parameters
        ----------
        coord : Coordinate
        label : string

        """
        if coord is _xtal.Site and label is None:
            self._pybind_value = None

        elif type(coord).__name__ is "Coordinate" or type(coord).__name__ is "MutableCoordinate":
            self._pybind_value = _xtal.Site(coord._pybind_value, label)

        else:
            self._pybind_value = _xtal.Site(coord, label)

    @classmethod
    def _from_pybind(cls, py_bind_value):
        """Returns a constructed _Site from
        a given _xtal.Site value

        Parameters
        ----------
        py_bind_value : _xtal.Site

        Returns
        -------
        _Site

        """
        value = cls(_xtal.Site, None)
        value._pybind_value = py_bind_value
        return value

    def cart(self):
        """Returns cartesian coordinates of the Site

        Returns
        -------
        np.array

        """
        return self._pybind_value._cart_const()

    def frac(self, lat):
        """Returns fractional coordinates of the Site

        Parameters
        ----------
        lat : Lattice

        Returns
        -------
        np.array

        """
        return self._pybind_value._frac_const(lat)

    def label(self):
        """Returns label of atom at the Site

        Returns
        -------
        string

        """
        return self._pybind_value._label_const()

    def set_compare_method(self, method, *args):
        """Determines what strategy to use for comparing
        Sites

        Parameters
        ----------
        method : Functor class that performs the comparision
        *args : Arguments needed to construct the Functor

        """
        self._equals = method(self, *args)

    def __eq__(self, other):
        """Passes the "other" to the cuurent compare functor
        and compares it to the self and returns the evaluation

        Parameters
        ----------
        other : Site or MutableSite

        Returns
        -------
        bool

        """
        if hasattr(self, '_equals') is False:
            self.set_compare_method(Equals, globaldef.tol)

        return self._equals(other)

    def __ne__(self, other):
        """Passes the "other" to the current compare functor
        and compares it to the self and returns opposite of the
        evaluation

        Parameters
        ----------
        other : Site or MutableSite

        Returns
        -------
        bool

        """
        return not self == other

    def __str__(self):
        """Returns the coordinate values along with the
        along with the atom type as a printable string

        Returns
        -------
        string

        """
        return self._pybind_value.__str__()

    def __rmul__(self, CartOp):
        """Applies the provided symmetry operation
        to the site and returns the transformed site

        Parameters
        ----------
        CartOp : cu.sym.CartOp

        Returns
        -------
        Site

        """
        return self._from_pybind(CartOp * self._pybind_value)


class Site(_Site):

    """Immutable Site Class. Defined as cartesian coordinates
    along with atom type. Handles all const site operations"""

    def __init__(self, coord, label):
        """
        Parameters
        ----------
        coord : Coordinate
        label : string

        """
        super().__init__(coord, label)


class MutableSite(_Site):

    """Mutable Site Class. Defined as cartesian coordinates along
    with atom type. Handles all non const site operations"""

    def __init__(self, coord, label):
        """
        Parameters
        ----------
        coord : MutableCoordinate
        label : string

        """
        super().__init__(coord, label)
