from __future__ import absolute_import

from ._xtal import *
from .single_block_wadsley_roth import *

def extra_function(self):
    print(self)
    print(make_niggli(self))
    return self.primitive()

Structure.extra_function=extra_function
