"""Place to collect miscellaneous useful functions."""

import numpy


def option_repr(value):
    """Return ``repr(value)`` safe for AstroWISP text option strings.

    The C++ tools are configured through ``boost::program_options`` using
    ``repr()`` of the Python values. NumPy 2 (NEP 51) reprs a numpy scalar as
    e.g. ``np.float64(1.0)`` instead of ``1.0``, which the option parser
    rejects. Converting numpy scalars to their native Python equivalent first
    yields a clean ``1.0`` while preserving the int/float distinction. Arrays
    and non-numpy values are returned to ``repr`` unchanged.
    """

    if isinstance(value, numpy.generic):
        value = value.item()
    return repr(value)


def flux_from_magnitude(magnitude, magnitude_1adu):
    """Return the flux corresponding to the given magnitude."""

    return 10.0**((magnitude_1adu - magnitude) / 2.5)
