from __future__ import absolute_import

from . import coordinate
from .coordinate import Coordinate
from .coordinate import MutableCoordinate
from .lattice import *
from .site import *
from .structure import *
from .globaldef import *
from ._xtal import make_niggli as _make_niggli
# from .single_block_wadsley_roth import *

def extra_function(self):
    print("I'm the extra function!")

def make_niggli(input_value):
    """Returns the niggli version of input_value. Type checks the argument
    to return a value with the same type.

    :input_value: casmutils.xtal.Lattice or casmutils.xtal.Structure
    :returns: casmutils.xtal.Lattice or casmutils.xtal.Structure 

    """
    if str(type(input_value))=="<class 'casmutils.xtal.structure.Structure'>": 
        return Structure._from_pybind(_make_niggli(input_value._pybind_value))
    elif str(type(input_value))=="<class 'casmutils.xtal.lattice.Lattice'>": 
        return Lattice._from_pybind(_make_niggli(input_value))
    else:
        return None

Coordinate.extra_function=extra_function
