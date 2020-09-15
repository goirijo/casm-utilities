from . import _mapping
from .structure import MappingReport
from ..xtal import Structure


def structure_score(mapping_data):
    """Calculates the individual lattice and basis scores from the ideal lattice,
    stretch tensor, and displacement matrix present in the provided MappingReport.
    Values are returned as a pair of [lattice_score, basis_score]

    Parameters
    ----------
    mapping_data : MappingReport

    Returns
    -------
    [float,float]

    """
    return _mapping.structure_score(mapping_data._pybind_value)


def map_structure(reference_struc, mapped_struc):
    """Using the default parameters for the StructureMapper, return a MappingReport report
    that maps one structure onto the other

    Parameters
    ----------
    reference_struc : xtal.Structure
    mapped_struc : xtal.Structure

    Returns
    -------
    MappingReport

    """
    return [
        MappingReport(r) for r in _mapping.map_structure(
            reference_struc._pybind_value, mapped_struc._pybind_value)
    ]
