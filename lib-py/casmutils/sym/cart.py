from __future__ import absolute_import

from . import _sym


class CartOp(_sym.CartOp):
    """Basic symmetry operation that deals with Cartesian representations.
    Conains a matrix, translation, and time reversal boolean."""

    # def __init__(self):
    #     """Initialize with the matrix, translation, and time reversal bool """
    #     _sym.CartOp.__init__(self)

    def __str__(self):
        """Concatenates __str__ of each member (matrix, vector, bool)
        Returns
        -------
        str

        """
        return self.matrix.__str__() + '\n\n' + self.translation.__str__(
        ) + '\n\n' + self.is_time_reversal_active.__str__()

    def __repr__(self):
        """Concatenates __repr__ of each member (matrix, vector, bool)
        Returns
        -------
        str

        """
        return self.matrix.__repr__() + '\n\n' + self.translation.__repr__(
        ) + '\n\n' + self.is_time_reversal_active.__repr__()

    def __iter__(self):
        yield from [self.matrix, self.translation, self.is_time_reversal_active]

    @classmethod
    def identity(cls):
        """Return an operation that has identity matrix, no translation, and no time reversal
        Returns
        -------
        CartOp

        """
        S=super().identity()
        return CartOp(S.matrix,S.translation,S.is_time_reversal_active)

    @classmethod
    def time_reversal(cls):
        """Return an operation that has identity matrix, no translation, and time reversal
        Returns
        -------
        CartOp

        """
        S=super().time_reversal()
        return CartOp(S.matrix,S.translation,S.is_time_reversal_active)

    @classmethod
    def translation_operation(cls, translation):
        """Return an operation that has identity matrix, the specified translation, and no time reversal

        Parameters
        ----------
        translation : 3x1 array

        Returns
        -------
        CartOp

        """
        S=super().translation_operation(translation)
        return CartOp(S.matrix,S.translation,S.is_time_reversal_active)

    @classmethod
    def point_operation(cls, matrix):
        """Return an operation that has the specified matrix, no translation, and no time reversal

        Parameters
        ----------
        matrix : 3x3 array

        Returns
        -------
        CartOp

        """
        S=super().point_operation(matrix)
        return CartOp(S.matrix,S.translation,S.is_time_reversal_active)
