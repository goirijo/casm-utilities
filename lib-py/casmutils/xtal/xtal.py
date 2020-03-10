from __future__ import absolute_import

# from ._xtal import *
from .coordinate import *
from .lattice import *
# from .single_block_wadsley_roth import *

def extra_function(self):
    print("I'm the extra function!")

Coordinate.extra_function=extra_function
