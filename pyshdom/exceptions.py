"""
Custom exceptions for pyshdom that are tailored to specific cases.
"""
class NegativeValueError(Exception):
    """Raised if a Value is identified as negative when it
    should not be (e.g. extinction). see. checks.check_positivity."""

class OutOfRangeError(Exception):
    """Raised if a Value is identified as being out of a required range
    when it should not be (e.g. single-scatter-albedo > 1.0).
    See checks.check_range."""

class MissingDimensionError(Exception):
    """Raised if a Variable in an xr.Dataset was expected
    to have a dimension with a specific name but it was not found.
    (e.g. the spatial coordinates should be labeled 'x', 'y', 'z').
    This ensures that priviledged names are present.
    See checks.check_hasdim"""

class LegendreTableError(Exception):
    """Raised if the formatting or Values in an xr.Dataset
    are not consistent with expectations for LegendreTable variables.
    See checks.check_legendre.
    """

class GridError(Exception):
    """Raised if an xr.Dataset has coordinates which are not consistent
    with the requirements of an SHDOM grid. See checks.check_grid.
    """

class SHDOMError(Exception):
    """Raised if there is an error that would call 'STOP' in SHDOM within
    the python wrapping, e.g. if there is an issue in memory allocation.
    """
