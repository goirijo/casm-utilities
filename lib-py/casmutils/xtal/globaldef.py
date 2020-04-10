#!@PYTHON@

""" All the global definitions are declared here """

"""Tolerance"""
tol = 1e-5

def _is_pybind_value(value):
    """Checks if the value passed
    is of _py_bind type

    Parameters
    ----------
    value : any type

    Returns
    -------
    bool

    """
    if ("_xtal" in type(value).__module__) is True:
        return True

    else:
        return False
