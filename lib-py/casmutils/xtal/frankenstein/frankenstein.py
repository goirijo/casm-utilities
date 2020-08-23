from __future__ import absolute_import

from ._frankenstein import translate_basis as _translate_basis
from .. import xtal

def translate_basis(struc,shift):
    """Create a new structure that has had each of its
    basis sites translated by the specified shift value.

    Parameters
    ----------
    struc : xtal.Structure
    shift : 1x3 np.array

    Returns
    -------
    TODO

    """
    return xtal.Structure._from_pybind(_translate_basis(struc._pybind_value,shift))
