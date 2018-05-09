"""Interface to the SuperPhot background library."""

from ctypes import\
    cdll,\
    c_double,\
    c_size_t,\
    c_void_p,\
    c_ulong,\
    c_bool,\
    c_char,\
    c_uint
from ctypes.util import find_library
import numpy.ctypeslib

#Naming convention imitated the one by ctypes.
#Type checking place holders require no content.
#pylint: disable=invalid-name
#pylint: disable=too-few-public-methods
class _c_background_extractor_p(c_void_p):
    pass

class _c_core_image_p(c_void_p):
    pass
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

def _initialize_library():
    """Prepare the SuperPhot background library for use."""

    library_fname = find_library('superphotbackground')
    if library_fname is None:
        raise OSError('Unable to find the SuperPhot background library.')
    library = cdll.LoadLibrary(library_fname)


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

    library.create_background_extractor.argtypes = [c_double,
                                                    c_double,
                                                    c_double,
                                                    _c_core_image_p,
                                                    c_double]
    library.create_background_extractor.restype = _c_background_extractor_p


    library.destroy_background_extractor.argtypes = [
        library.create_background_extractor.restype
    ]
    library.destroy_background_extractor.restype = None


    library.add_source_list_to_background_extractor.argtypes = [
        library.create_background_extractor.restype,
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=c_double,
                                  ndim=1,
                                  flags='C_CONTIGUOUS'),
        c_size_t
    ]
    library.add_source_list_to_background_extractor.restype = None

    library.get_all_backgrounds.argtypes = [
        library.create_background_extractor.restype,
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_double,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
        ndpointer_or_null(dtype=c_uint,
                          ndim=1,
                          flags='C_CONTIGUOUS'),
    ]
    library.get_all_backgrounds.restype = None

    return library

class BackgroundExtractor:
    """
    Measure the background level for each source in an image.

    Attributes:
        image:    The image being processed.

        inner_radius:    The size of the aperture aronud each source within
            which pixels are excluded from background measurement.

        outer_radius:    The outer rim of the aperture around each source within
            which unrejected pixels are included in the background measurement.

        error_confidence:    The confidence level to use for estimating the
            background error.
    """

    _library = _initialize_library()

    def __init__(self,
                 image,
                 inner_radius,
                 outer_radius,
                 error_confidence=0.68):
        """
        Create a background extractor with the given parameters.

        Args: see class attributes.

        Returns: None
        """

        self.image = image
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.error_confidence = error_confidence
        self._library_image = self._library.create_core_image(
            image.shape[1],
            image.shape[0],
            image,
            None,
            None,
            True
        )
        self._library_extractor = self._library.create_background_extractor(
            inner_radius,
            outer_radius,
            inner_radius,
            self._library_image,
            error_confidence
        )

    def __call__(self, source_x, source_y):
        """
        Measure the background under the sources with the given coordinates.

        Args:
            source_x:    The `x` coordinates of the sources within the image.

            source_y:    The `y` coordinates of the sources within the image.

        Returns:
            tuple:
                numpy.array:
                    The estimate of the background under each source in the same
                    order as the input sources.

                numpy.array:
                    The estimate of the uncertainty in the background under each
                    source in the same order as the input sources.

                numpy.array:
                    The number of pixels which were used to derive the
                    background and its uncertainty.
        """

        assert source_x.size == source_y.size

        self._library.add_source_list_to_background_extractor(
            self._library_extractor,
            source_x,
            source_y,
            source_x.size
        )

        bg_value = numpy.empty(source_x.size, dtype=c_double)
        bg_error = numpy.empty(source_x.size, dtype=c_double)
        bg_numpix = numpy.empty(source_x.size, dtype=c_uint)
        self._library.get_all_backgrounds(
            self._library_extractor,
            bg_value,
            bg_error,
            bg_numpix
        )
        return bg_value, bg_error, bg_numpix

    def __del__(self):
        r"""Destroy the image and extractor created in :meth:`__init__`\ ."""

        self._library.destroy_core_image(self._library_image)
        self._library.destroy_background_extractor(self._library_extractor)
