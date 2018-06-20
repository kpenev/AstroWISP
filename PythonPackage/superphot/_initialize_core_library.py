"""Define function preparing SuperPhot Core functions for use."""

from ctypes import c_ulong, c_char, c_bool, c_double, c_void_p
import numpy.ctypeslib

#Naming convention imitates the one by ctypes.
#pylint: disable=invalid-name

#Type checking place holders require no content.
#pylint: disable=too-few-public-methods
class _c_core_image_p(c_void_p):
    """Placeholder for CoreImage opaque struct."""

class _c_core_sub_pixel_map_p(c_void_p):
    """Placeholder for CoreSubPixelMap opaque struct."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def ndpointer_or_null(*args, **kwargs):
    """
    Allow None (->NULL) to be passed for c-style array function arguments.

    Modified from:
    http://stackoverflow.com/questions/32120178/how-can-i-pass-null-to-an-external-library-using-ctypes-with-an-argument-decla
    """

    base = numpy.ctypeslib.ndpointer(*args, **kwargs)

    #Call signature dictated by numpy.ctypeslib
    #pylint: disable=unused-argument
    def from_param(cls, obj):
        """Construct numpy.ndpointer from the given object."""

        if obj is None:
            return obj
        return base.from_param(obj)
    #pylint: enable=unused-argument

    return type(base.__name__,
                (base,),
                {'from_param': classmethod(from_param)})

def _initialize_core_library(library):
    """Set-up the argument and return types of Core library functions."""

    library.create_core_image.argtypes = [
        c_ulong,
        c_ulong,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=2,
                                  flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=2,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_char,
                          ndim=2,
                          flags='C_CONTIGUOUS'),
        c_bool
    ]
    library.create_core_image.restype = _c_core_image_p

    library.destroy_core_image.argtypes = [
        library.create_core_image.restype
    ]
    library.destroy_core_image.restype = None

    library.create_core_subpixel_map.argtypes = [
        c_ulong,
        c_ulong,
        ndpointer_or_null(dtype=c_double,
                          ndim=2,
                          flags='C_CONTIGUOUS')
    ]
    library.create_core_subpixel_map.restype = _c_core_sub_pixel_map_p

    library.destroy_core_subpixel_map.argtypes = [
        library.create_core_subpixel_map.restype
    ]
    library.destroy_core_subpixel_map.restype = None
