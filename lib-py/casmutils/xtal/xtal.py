from __future__ import absolute_import

from . import coordinate
from .coordinate import Coordinate
from .coordinate import MutableCoordinate
from .lattice import *
from .site import *
from .structure import *
from .symmetry import *
from .globaldef import *
from ._xtal import make_niggli as _make_niggli
from ._xtal import make_superstructure as _make_superstructure
from ._xtal import make_primitive as _make_primitive
# from .single_block_wadsley_roth import *


def extra_function(self):
    print("I'm the extra function!")


def make_niggli(input_value):
    """Returns the niggli version of input_value. Type checks the argument
    to return a value with the same type.

    :input_value: casmutils.xtal.lattice.Lattice or casmutils.xtal.structure.Structure
    :returns: casmutils.xtal.lattice.Lattice or casmutils.xtal.structure.Structure 

    """
    if isinstance(input_value, Structure):
        return Structure._from_pybind(_make_niggli(input_value._pybind_value))
    elif isinstance(input_value, Lattice):
        return Lattice._from_pybind(_make_niggli(input_value))
    else:
        raise ValueError


def make_superstructure(structure, transformation_matrix):
    """Returns the superstructure of the given structure,
    scaling the lattice by the given transformation matrix

    :structure: casmutils.xtal.structure.Structure
    :transformation_matrix: np.array(int32[3,3])
    :returns: casmutils.xtal.structure.Structure

    """
    return Structure._from_pybind(
        _make_superstructure(structure._pybind_value, transformation_matrix))


def make_primitive(structure):
    """Returns the primitive version of the given structure.

    :structure: casmutils.xtal.structure.Structure
    :returns: casmutils.xtal.structure.Structure

    """
    return Structure._from_pybind(_make_primitive(structure._pybind_value))


Coordinate.extra_function = extra_function
