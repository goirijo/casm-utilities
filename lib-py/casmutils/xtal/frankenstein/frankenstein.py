from __future__ import absolute_import

from ._frankenstein import translate_basis as _translate_basis
from ._frankenstein import stack as _stack
from .. import xtal


def translate_basis(struc, shift):
    """Create a new structure that has had each of its
    basis sites translated by the specified shift value.

    Parameters
    ----------
    struc : xtal.Structure
    shift : 1x3 np.array

    Returns
    -------
    xtal.Structure

    """
    return xtal.Structure._from_pybind(
        _translate_basis(struc._pybind_value, shift))


def stack(strucs):
    """Given a series of structures that have the same ab values, stack them
    in the given order along the c direction to create a tower. If the ab values of
    all the structures aren't related by a rigid rotation, the structures will be
    distorted to match the first given structure.

    Parameters
    ----------
    strucs : xtal.Structure values in the stacking order, first argument ends up on the bottom

    Returns
    -------
    xtal.Structure

    """
    return xtal.Structure._from_pybind(
        _stack([s._pybind_value for s in strucs]))
