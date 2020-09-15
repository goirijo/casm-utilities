from . import _xtal
from ..sym.cart import CartOp
from .lattice import *
from .structure import *


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
    return [CartOp(op) for op in _xtal.make_point_group(lattice, tol)]


def make_factor_group(structure, tol):
    """Calculate all symmetry operations that map the given
    structure onto itself

    Parameters
    ----------
    structure : Structure
    tol : float

    Returns
    -------
    list(cu.sym.CartOp)

    """
    return [CartOp(op) for op in _xtal.make_factor_group(structure._pybind_value, tol)]


def symmetrize(lattice_or_structure, enforced_group):
    """Gives a symmetrized version of the input lattice/structure such 
    that it obeys the given enforced symmetry group

    :lattice_or_structure: Lattice or Structure
    :enforced_group: list(cu.sym.CartOp)
    :returns: Lattice or Structure

    """
    if isinstance(lattice_or_structure, Lattice):
        return _xtal._symmetrize_lattice(lattice_or_structure, enforced_group)
    elif isinstance(lattice_or_structure, Structure):
        return Structure._from_pybind(_xtal._symmetrize_structure(lattice_or_structure._pybind_value, enforced_group))
    else:
        raise ValueError("symmetrize only works on Structure or Lattice types")
