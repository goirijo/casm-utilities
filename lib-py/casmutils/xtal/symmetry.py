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
