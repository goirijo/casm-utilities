from __future__ import absolute_import

from ._xtal import *
from .single_block_wadsley_roth import *

def test_function(arg1, arg2):
    """This function does nothing, I'm just testing sphinx.

    Parameters
    ----------
    arg1 : str
    arg2 : str

    Returns
    -------
    Null

    """
    return

def extra_function(self):
    print(self)
    print(make_niggli(self))
    return self.primitive()

Structure.extra_function=extra_function
