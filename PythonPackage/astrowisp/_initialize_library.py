"""Load the C library and define the function interface."""

from os import path

from ctypes import (
    cdll,
    c_double,
    c_size_t,
    c_void_p,
    c_ulong,
    c_bool,
    c_char,
    c_uint,
    c_char_p,
    c_long,
    c_byte,
    POINTER,
)
from ctypes.util import find_library
import numpy.ctypeslib

# Naming convention imitates the one by ctypes.
# pylint: disable=invalid-name


# Type checking place holders require no content.
# pylint: disable=too-few-public-methods
class _c_core_image_p(c_void_p):
    """Placeholder for CoreImage opaque struct."""


class _c_core_sub_pixel_map_p(c_void_p):
    """Placeholder for CoreSubPixelMap opaque struct."""


class _c_background_extractor_p(c_void_p):
    """Placeholder for BackgroundExtractor opaque struct."""


class _c_h5io_data_tree_p(c_void_p):
    """Placeholder for pointer to H5IODataTree opaque struct."""


class _c_fitting_configuration(c_void_p):
    """Placeholder for the FittingConfiguration opaque struct."""


class _c_subpixphot_configuration(c_void_p):
    """Placeholder for the SubPixPhotConfiguration opaque struct."""


class _c_piecewise_bicubic_psf_map_p(c_void_p):
    """Placeholder for the PiecewiseBicubicPSFMap opaque struct."""


class _c_piecewise_bicubic_psf_p(c_void_p):
    """Placeholder for the PiecewiseBicubicPSF opaque struct."""


# pylint: enable=invalid-name
# pylint: enable=too-few-public-methods


def ndpointer_or_null(*args, **kwargs):
    """
    Allow None (->NULL) to be passed for c-style array function arguments.

    Modified from:
    http://stackoverflow.com/questions/32120178/how-can-i-pass-null-to-an-external-library-using-ctypes-with-an-argument-decla
    """

    base = numpy.ctypeslib.ndpointer(*args, **kwargs)

    # Call signature dictated by numpy.ctypeslib
    # pylint: disable=unused-argument
    def from_param(cls, obj):
        """Construct numpy.ndpointer from the given object."""

        if obj is None:
            return obj
        return base.from_param(obj)

    # pylint: enable=unused-argument

    return type(base.__name__, (base,), {"from_param": classmethod(from_param)})


def _setup_core_interface(library):
    """Set-up the argument and return types of Core library functions."""

    library.create_core_image.argtypes = [
        c_ulong,
        c_ulong,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
        ndpointer_or_null(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
        ndpointer_or_null(dtype=c_char, ndim=2, flags="C_CONTIGUOUS"),
        c_bool,
    ]
    library.create_core_image.restype = _c_core_image_p

    library.destroy_core_image.argtypes = [library.create_core_image.restype]
    library.destroy_core_image.restype = None

    library.create_core_subpixel_map.argtypes = [
        c_ulong,
        c_ulong,
        ndpointer_or_null(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
    ]
    library.create_core_subpixel_map.restype = _c_core_sub_pixel_map_p

    library.destroy_core_subpixel_map.argtypes = [
        library.create_core_subpixel_map.restype
    ]
    library.destroy_core_subpixel_map.restype = None


def _setup_io_interface(library):
    """Set-up the argument and return types of the I/O library functions."""

    library.create_result_tree.argtypes = [c_void_p, c_char_p]
    library.create_result_tree.restype = _c_h5io_data_tree_p

    library.destroy_result_tree.argtypes = [library.create_result_tree.restype]
    library.destroy_result_tree.restype = None

    library.query_result_tree.argtypes = [
        library.create_result_tree.restype,
        c_char_p,
        c_char_p,
        c_void_p,
    ]
    library.query_result_tree.restype = c_bool

    library.get_psf_map_variables.argtypes = [
        library.create_result_tree.restype,
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
    ]
    library.get_psf_map_variables.restype = None

    library.list_tree_quantities.argtypes = [
        library.create_result_tree.restype,
        POINTER(POINTER(c_char_p)),
    ]
    library.list_tree_quantities.restype = c_uint

    library.parse_hat_mask.argtypes = [
        c_char_p,
        c_long,
        c_long,
        numpy.ctypeslib.ndpointer(dtype=c_byte, ndim=1, flags="C_CONTIGUOUS"),
    ]
    library.parse_hat_mask.restype = None

    library.mask_flags = {
        flag: c_byte.in_dll(library, "MASK_" + flag).value
        for flag in [
            "OK",
            "CLEAR",
            "FAULT",
            "HOT",
            "COSMIC",
            "OUTER",
            "OVERSATURATED",
            "LEAKED",
            "SATURATED",
            "INTERPOLATED",
            "BAD",
            "ALL",
            "NAN",
        ]
    }

    library.update_result_tree.argtypes = [
        # Which quantity to add/overwrite.
        c_char_p,
        # The value(s) to set for the quantity.
        c_void_p,
        # Data type of quantity.
        c_char_p,
        # The number of entries in the quantity
        c_uint,
        # The tree to update
        library.create_result_tree.restype,
    ]
    library.update_result_tree.restype = None

    library.set_psf_map_variables.argtypes = [
        # variblae names
        POINTER(c_char_p),
        # variable values
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
        # number of variables
        c_uint,
        # number of sources
        c_uint,
        # The index of the image to set the PSF map for.
        c_uint,
        # The tree to update
        library.create_result_tree.restype,
    ]
    library.set_psf_map_variables.restype = None

    library.export_free.argtypes = [c_void_p]
    library.export_free.restype = None


def _setup_background_interface(library):
    """Set-up the argument and return types of the background library funcs."""

    library.create_background_extractor.argtypes = [
        c_double,
        c_double,
        c_double,
        library.create_core_image.restype,
        c_double,
    ]
    library.create_background_extractor.restype = _c_background_extractor_p

    library.destroy_background_extractor.argtypes = [
        library.create_background_extractor.restype
    ]
    library.destroy_background_extractor.restype = None

    library.add_source_list_to_background_extractor.argtypes = [
        library.create_background_extractor.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        c_size_t,
    ]
    library.add_source_list_to_background_extractor.restype = None

    library.get_all_backgrounds.argtypes = [
        library.create_background_extractor.restype,
        ndpointer_or_null(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        ndpointer_or_null(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        ndpointer_or_null(dtype=c_uint, ndim=1, flags="C_CONTIGUOUS"),
    ]
    library.get_all_backgrounds.restype = None

    return library


def _setup_psf_interface(library):
    """Set-up the argument and return types of the PSF library funcs."""

    library.create_piecewise_bicubic_psf_map.argtypes = [
        library.create_result_tree.restype
    ]
    library.create_piecewise_bicubic_psf_map.restype = _c_piecewise_bicubic_psf_map_p

    library.destroy_piecewise_bicubic_psf_map.argtypes = [
        library.create_piecewise_bicubic_psf_map.restype
    ]
    library.destroy_piecewise_bicubic_psf_map.restype = None

    library.evaluate_piecewise_bicubic_psf_map.argtypes = [
        library.create_piecewise_bicubic_psf_map.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
    ]
    library.evaluate_piecewise_bicubic_psf_map.restype = _c_piecewise_bicubic_psf_p

    library.destroy_piecewise_bicubic_psf.argtypes = [
        library.evaluate_piecewise_bicubic_psf_map.restype
    ]
    library.destroy_piecewise_bicubic_psf.restype = None

    library.evaluate_piecewise_bicubic_psf.argtypes = [
        library.evaluate_piecewise_bicubic_psf_map.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
    ]
    library.evaluate_piecewise_bicubic_psf.restype = None

    library.integrate_piecewise_bicubic_psf.argtypes = [
        library.evaluate_piecewise_bicubic_psf_map.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
        c_uint,
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags="C_CONTIGUOUS"),
    ]
    library.integrate_piecewise_bicubic_psf.restype = None


def _setup_fitpsf_interface(library):
    """Set-up the argument and return types of the PSF fitting library funcs."""

    library.create_psffit_configuration.argtypes = []
    library.create_psffit_configuration.restype = _c_fitting_configuration

    library.destroy_psffit_configuration.argtypes = [
        library.create_psffit_configuration.restype
    ]

    library.update_psffit_configuration.argtypes = [
        c_bool,
        library.create_psffit_configuration.restype,
    ]
    library.update_psffit_configuration.restype = None

    library.piecewise_bicubic_fit.argtypes = [
        # pixel_values
        POINTER(POINTER(c_double)),
        # pixel_errors
        POINTER(POINTER(c_double)),
        # pixel_masks
        POINTER(POINTER(c_char)),
        # number_images
        c_ulong,
        # image_x_resolution
        c_ulong,
        # image_y_resolution
        c_ulong,
        # source_ids
        POINTER(POINTER(c_char_p)),
        # source_coordinates
        POINTER(POINTER(c_double)),
        # psf_terms
        POINTER(POINTER(c_double)),
        # enabled
        POINTER(POINTER(c_bool)),
        # number_sources
        numpy.ctypeslib.ndpointer(dtype=c_ulong, ndim=1, flags="C_CONTIGUOUS"),
        # number_terms
        c_ulong,
        # backgrounds
        library.create_background_extractor.restype,
        # configuration
        library.create_psffit_configuration.restype,
        # subpixel_sensitivities
        numpy.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags="C_CONTIGUOUS"),
        # subpix_x_resolution
        c_ulong,
        # subpix_y_resolution
        c_ulong,
        # ouput_data_tree
        library.create_result_tree.restype,
    ]
    library.piecewise_bicubic_fit.restype = c_bool


def _setup_subpixphot_interface(library):
    """Set-up the argument and return tyes of the aperture photometry funcs."""

    library.create_subpixphot_configuration.argtypes = []
    library.create_subpixphot_configuration.restype = _c_subpixphot_configuration

    library.destroy_subpixphot_configuration.argtypes = [
        library.create_subpixphot_configuration.restype
    ]
    library.destroy_subpixphot_configuration.restype = None

    library.update_subpixphot_configuration.argtypes = [
        library.create_subpixphot_configuration.restype
    ]
    library.update_subpixphot_configuration.restype = None

    library.subpixphot.argtypes = [
        _c_core_image_p,
        _c_core_sub_pixel_map_p,
        library.create_subpixphot_configuration.restype,
        library.create_result_tree.restype,
        c_uint,
    ]
    library.subpixphot.restype = None


def _initialize_library():
    """Prepare the astrowisp library for use."""

    lib_fname = path.join(path.dirname(__file__), "libastrowisp.")
    library = None
    for ext in ["so", "dylib", "dll"]:
        if path.exists(lib_fname + ext):
            assert library is None
            library = cdll.LoadLibrary(lib_fname + ext)
    if library is None:
        lib_fname = find_library("astrowisp")
        assert lib_fname is not None
        library = cdll.LoadLibrary(lib_fname)

    _setup_core_interface(library)
    _setup_io_interface(library)
    _setup_background_interface(library)
    _setup_psf_interface(library)
    _setup_fitpsf_interface(library)
    _setup_subpixphot_interface(library)

    return library


def get_astrowisp_library():
    """Return the shared astrowisp library."""

    if not hasattr(get_astrowisp_library, "result"):
        get_astrowisp_library.result = _initialize_library()

    return get_astrowisp_library.result
