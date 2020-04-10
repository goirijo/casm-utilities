from . import _mapping

MappingReport = _mapping.MappingReport

class MappingInput(_mapping.MappingInput):
    """Specifies all possible mapping specifications except for the reference
    structure itself, and the symmetry operations that should be considered."""

    def __init__(self,
                 strain_weight=None,
                 max_volume_change=None,
                 options=None,
                 tol=None,
                 min_vacancy_fraction=None,
                 max_vacancy_fraction=None,
                 k_best_maps=None,
                 max_cost=None,
                 min_cost=None,
                 keep_invalid_mapping_nodes=None,
                 impose_reference_lattice=None,
                 assume_ideal_structure=None,
                 assume_ideal_lattice=None):
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
            self.strain_weight = strain_weight

        if max_volume_change is not None:
            self.max_volume_change = max_volume_change

        #TODO: Implement as string
        if options is not None:
            raise NotImplemented(
                "Specifiying options as string not yet implemented,\
                    you can override the value manually post construction with the appropriate\
                    integer value post construction (sorry)")
            self.options = options

        if tol is not None:
            self.tol = tol

        if min_vacancy_fraction is not None:
            self.min_vacancy_fraction = min_vacancy_fraction

        if max_vacancy_fraction is not None:
            self.max_vacancy_fraction = max_vacancy_fraction

        if k_best_maps is not None:
            self.k_best_maps = k_best_maps

        if max_cost is not None:
            self.max_cost = max_cost

        if min_cost is not None:
            self.min_cost = min_cost

        if keep_invalid_mapping_nodes is not None:
            self.keep_invalid_mapping_nodes = keep_invalid_mapping_nodes

        if impose_reference_lattice is not None:
            self.impose_reference_lattice = impose_reference_lattice

        if assume_ideal_structure is not None:
            self.assume_ideal_structure = assume_ideal_structure

        if assume_ideal_lattice is not None:
            self.assume_ideal_lattice = assume_ideal_lattice

    def __str__(self):
        as_str=""

        as_str+="strain_weight:\n"+str(self.strain_weight)+"\n\n"
        as_str+="max_volume_change:\n"+str(self.max_volume_change)+"\n\n"
        as_str+="options:\n"+str(self.options)+"\n\n"
        as_str+="tol:\n"+str(self.tol)+"\n\n"
        as_str+="min_vacancy_fraction:\n"+str(self.min_vacancy_fraction)+"\n\n"
        as_str+="max_vacancy_fraction:\n"+str(self.max_vacancy_fraction)+"\n\n"
        as_str+="k_best_maps:\n"+str(self.k_best_maps)+"\n\n"
        as_str+="max_cost:\n"+str(self.max_cost)+"\n\n"
        as_str+="min_cost:\n"+str(self.min_cost)+"\n\n"
        as_str+="keep_invalid_mapping_nodes:\n"+str(self.keep_invalid_mapping_nodes)+"\n\n"
        as_str+="impose_reference_lattice:\n"+str(self.impose_reference_lattice)+"\n\n"
        as_str+="assume_ideal_structure:\n"+str(self.assume_ideal_structure)+"\n\n"
        as_str+="assume_ideal_lattice:\n"+str(self.assume_ideal_lattice)

        return as_str

class StructureMapper(_mapping.StructureMapper_f):

    """Once constructed with a reference structure and the relevant mapping parameters,
    use the call operator to map other structures onto the initial reference structure.
    See the MappingInput class for the possible parameters that can be set when mapping
    structures. 
    
    TODO: Add more info"""

    def __init__(self,reference_structure, mapping_input=None, point_group=[], allowed_species=[], **kwargs):
        """The reference structure is always required. If no other values are given,
        then default values will be created for the MappingInput. A custom symmetry
        group (which would usually be the point group of the reference structure) can
        be sepecified. When mapping structures to a reference site that should allow
        multiple species types, the allowed_species parameter can specify which types
        are allowed to occupy each site in the reference structure.

        Parameters
        ----------
        reference_structure : xtal.Structure
        mapping_input : MappingInput
        point_group : List[sym.CartOp]
        allowed_species : List[List[str]]
        kwargs : values to specify for the MappingInput, ignored if mapping_input is provided


        """
        self._reference_structure = reference_structure

        #TODO: Ideally you can pass both MappingInput and kwargs, which override field by field?
        if mapping_input is None:

            if "k" in kwargs:
                kwargs["k_best_maps"]=kwargs["k"]
                del kwargs["k"]

            mapping_input=MappingInput(**kwargs)

        self._mapping_input=mapping_input
        super().__init__(reference_structure, mapping_input, point_group, allowed_species)

