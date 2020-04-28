from . import _xtal
from ..sym.cart import CartOp

def make_point_group(lattice, tol):
    """Calculate all symmetry operations that map the given
    lattice onto itself

    Parameters
    ----------
    lattice : Lattice
    tol : float

    Returns
    -------
    list(cu.sym.CartOp)

    """
    return [CartOp(op) for op in _xtal.make_point_group(lattice,tol)]

def make_factor_group(structure, tol):
    """Calculate all symmetry operations that map the given
    structure onto itself

    Parameters
    ----------
    lattice : Structure
    tol : float

    Returns
    -------
    list(cu.sym.CartOp)

    """
    return [CartOp(op) for op in _xtal.make_factor_group(structure._pybind_value,tol)]

def symmetrize(lattice, point_group):
    """Gives a symmetrized version of the input lattice such 
    that it obeys the given point group

    :lattice: Lattice
    :point_group: list(cu.sym.CartOp)
    :returns: Lattice

    """
    return _xtal.symmetrize(lattice, point_group)
