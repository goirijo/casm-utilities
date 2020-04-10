from __future__ import absolute_import

# from ._xtal import *
from . import coordinate
from .coordinate import Coordinate
from .coordinate import MutableCoordinate
from .lattice import *
from .site import *
from .structure import *
from .globaldef import *

# from .single_block_wadsley_roth import *

def extra_function(self):
    print("I'm the extra function!")

Coordinate.extra_function=extra_function
