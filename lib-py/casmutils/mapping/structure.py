from . import _mapping


class MappingReport():
    """Describes the results of mapping one structure onto another, including
    a description of the tensors required to match the lattices, permutations
    and displacements required to map the basis, and both individual costs associated
    with each part of the mapping, as well as a combined lattice and basis mapping score"""

    def __init__(self, pybind_value):
        self._pybind_value = pybind_value

        self.isometry = self._pybind_value.isometry
        self.stretch = self._pybind_value.stretch
        self.translation = self._pybind_value.translation
        self.displacement = self._pybind_value.displacement
        self.permutation = self._pybind_value.permutation
        self.lattice_cost = self._pybind_value.lattice_cost
        self.basis_cost = self._pybind_value.basis_cost
        self.cost = self._pybind_value.cost
        self.reference_lattice = self._pybind_value.reference_lattice
        self.mapped_lattice = self._pybind_value.mapped_lattice

    def __str__(self):
        as_str = ""

        as_str += "isometry:\n"
        as_str += self.isometry.__str__()
        as_str += "\n\n"

        as_str += "stretch:\n"
        as_str += self.stretch.__str__()
        as_str += "\n\n"

        as_str += "translation:\n"
        as_str += self.translation.__str__()
        as_str += "\n\n"

        as_str += "displacement:\n"
        as_str += self.displacement.__str__()
        as_str += "\n\n"

        as_str += "permutation:\n"
        as_str += self.permutation.__str__()
        as_str += "\n\n"

        as_str += "lattice cost:\n"
        as_str += str(self.lattice_cost)
        as_str += "\n\n"

        as_str += "basis cost:\n"
        as_str += str(self.basis_cost)
        as_str += "\n\n"

        as_str += "cost:\n"
        as_str += str(self.cost)
        as_str += "\n\n"

        as_str += "reference lattice:\n"
        as_str += self.reference_lattice.__str__()
        as_str += "\n\n"

        as_str += "mapped lattice:\n"
        as_str += self.mapped_lattice.__str__()
        as_str += "\n\n"

        return as_str


class MappingInput(_mapping.MappingInput):
    """Specifies all possible mapping specifications except for the reference
    structure itself, and the symmetry operations that should be considered."""

    @classmethod
    def _sanitize_kwargs(cls, kwargs):
        """Some options will have more than one way to be specified,
        this function will translate shorthands to the most explicit
        keyword available.

        Parameters
        ----------
        kwargs : dict or MappingInput input parameters

        Returns
        -------
        dict

        """
        if "k" in kwargs:
            kwargs["k_best_maps"] = kwargs["k"]
            del kwargs["k"]

        #TODO: You could potentially set defaults here too? Probably
        #better to leave that in the c++ implementatiton though
        return kwargs

    def __setattr__(self,name,value):
        if not hasattr(self,name):
            raise AttributeError("{} is not a valid keyword for MappingInput".format(name))
        super().__setattr__(name,value)

    def __init__(self, **kwargs):
        #TODO: list what the default values are here?

        """Specify only values for which you don't want to keep the defaults.

        Parameters
        ----------
        strain_weight : float, optional
        max_volume_change : float, optional
        options : str, optional
        tol : float, optional
        min_vacancy_fraction : float, optional
        max_vacancy_fraction : float, optional
        k_best_maps : int, optional
        max_cost : float, optional
        min_cost : float, optional
        keep_invalid_mapping_nodes : bool, optional
        impose_reference_lattice : bool, optional
        assume_ideal_structure : bool, optional
        assume_ideal_lattice : bool, optional
        use_crystal_symmetry : bool, optional

        """
        _mapping.MappingInput.__init__(self)
        kwargs=MappingInput._sanitize_kwargs(kwargs)

        for k in kwargs:
            self.__setattr__(k,kwargs[k])

    def __str__(self):
        as_str = ""

        as_str += "strain_weight:\n" + str(self.strain_weight) + "\n\n"
        as_str += "max_volume_change:\n" + str(self.max_volume_change) + "\n\n"
        as_str += "options:\n" + str(self.options) + "\n\n"
        as_str += "tol:\n" + str(self.tol) + "\n\n"
        as_str += "min_vacancy_fraction:\n" + str(
            self.min_vacancy_fraction) + "\n\n"
        as_str += "max_vacancy_fraction:\n" + str(
            self.max_vacancy_fraction) + "\n\n"
        as_str += "k_best_maps:\n" + str(self.k_best_maps) + "\n\n"
        as_str += "max_cost:\n" + str(self.max_cost) + "\n\n"
        as_str += "min_cost:\n" + str(self.min_cost) + "\n\n"
        as_str += "keep_invalid_mapping_nodes:\n" + str(
            self.keep_invalid_mapping_nodes) + "\n\n"
        as_str += "impose_reference_lattice:\n" + str(
            self.impose_reference_lattice) + "\n\n"
        as_str += "assume_ideal_structure:\n" + str(
            self.assume_ideal_structure) + "\n\n"
        as_str += "assume_ideal_lattice:\n" + str(
            self.assume_ideal_lattice) + "\n\n"
        as_str += "use_crystal_symmetry:\n" + str(self.use_crystal_symmetry)

        return as_str


class StructureMapper:
    """Once constructed with a reference structure and the relevant mapping parameters,
    use the call operator to map other structures onto the initial reference structure.
    See the MappingInput class for the possible parameters that can be set when mapping
    structures.
    """

    # TODO: Add more info

    def __init__(self,
                 reference_structure,
                 mapping_input=None,
                 factor_group=[],
                 allowed_species=[],
                 **kwargs):
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
        factor_group : List[sym.CartOp]
        allowed_species : List[List[str]]
        kwargs : values to specify for the MappingInput, ignored if mapping_input is provided


        """
        self._reference_structure = reference_structure

        #TODO: Ideally you can pass both MappingInput and kwargs, which override field by field?
        if mapping_input is None:

            mapping_input = MappingInput(**kwargs)

        self._mapping_input = mapping_input
        self._pybind_value = _mapping.StructureMapper_f(
            reference_structure._pybind_value, mapping_input, factor_group,
            allowed_species)

    def __call__(self, structure):
        return [
            MappingReport(r)
            for r in self._pybind_value(structure._pybind_value)
        ]
