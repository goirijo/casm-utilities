from . import _mapping

MappingReport=_mapping.MappingReport

class MappingInput(_mapping.MappingInput):

    """Specifies all possible mapping specifications except for the reference
    structure itself, and the symmetry operations that should be considered."""

    def __init__(self, strain_weight=None, max_volume_change=None, options=None, tol=None, min_va_fraction=None, max_va_fraction=None, k_best_maps=None, max_cost=None, min_cost=None, keep_invalid_mapping_nodes=None, impose_reference_lattice=None, assume_ideal_structure=None, assume_ideal_lattice=None):
        """Specify only values for which you don't want to keep the defaults.

        Parameters
        ----------
        strain_weight : float, optional
        max_volume_change : float, optional
        options : str, optional
        tol : float, optional
        min_va_fraction : float, optional
        max_va_fraction : float, optional
        k_best_maps : int, optional
        max_cost : float, optional
        min_cost : float, optional
        keep_invalid_mapping_nodes : bool, optional
        impose_reference_lattice : bool, optional
        assume_ideal_structure : bool, optional
        assume_ideal_lattice : bool, optional

        """
        _mapping.MappingInput.__init__(self)

        if strain_weight is not None:
            self.strain_weight=strain_weight

        if max_volume_change is not None:
            self.max_volume_change=max_volume_change

        #TODO: Implement as string
        if options is not None:
            raise NotImplemented("Specifiying options as string not yet implemented,\
                    you can override the value manually post construction with the appropriate\
                    integer value post construction (sorry)")
            self.options=options

        if tol is not None:
            self.tol=tol

        if min_va_fraction is not None:
            self.min_va_fraction=min_va_fraction

        if max_va_fraction is not None:
            self.max_va_fraction=max_va_fraction

        if k_best_maps is not None:
            self.k_best_maps=k_best_maps

        if max_cost is not None:
            self.max_cost=max_cost

        if min_cost is not None:
            self.min_cost=min_cost

        if keep_invalid_mapping_nodes is not None:
            self.keep_invalid_mapping_nodes=keep_invalid_mapping_nodes

        if impose_reference_lattice is not None:
            self.impose_reference_lattice=impose_reference_lattice

        if assume_ideal_structure is not None:
            self.assume_ideal_structure=assume_ideal_structure

        if assume_ideal_lattice is not None:
            self.assume_ideal_lattice=assume_ideal_lattice

