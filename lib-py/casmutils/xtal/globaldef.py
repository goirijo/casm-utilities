#!@PYTHON@

""" All the global definitions are declared here """

"""Tolerance"""
tol = 1e-5

def _is_pybind_value(value):
    """Checks if the value passed
    is of py_bind type

    Parameters
    ----------
    value : any type

    Returns
    -------
    bool

    """
    if type(value).__name__ =="Site" and hasattr(value,"_cart_const"):
        return True

    elif type(value).__name__=="Coordinate" and hasattr(value,"_cart_const"):
        return True

    else:
        return False
