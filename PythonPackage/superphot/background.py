"""Interface to the SuperPhot background library."""

from ctypes import\
    cdll,\
    c_double,\
    c_size_t,\
    c_void_p,\
    c_uint
from ctypes.util import find_library
import numpy.ctypeslib

from superphot._initialize_core_library import\
    _initialize_core_library,\
    ndpointer_or_null

#Naming convention imitates the one by ctypes.
#pylint: disable=invalid-name

#Type checking place holders require no content.
#pylint: disable=too-few-public-methods
class _c_background_extractor_p(c_void_p):
    """Placeholder for BackgroundExtractor opaque struct."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def _initialize_library():
    """Prepare the SuperPhot background library for use."""

    library_fname = find_library('superphotbackground')
    if library_fname is None:
        raise OSError('Unable to find the SuperPhot background library.')
    library = cdll.LoadLibrary(library_fname)

    _initialize_core_library(library)

    library.create_background_extractor.argtypes = [
        c_double,
        c_double,
        c_double,
        library.create_core_image.restype,
        c_double
    ]
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

#The __init__, __del__ and __call__ methods justify making this a class.
#pylint: disable=too-few-public-methods
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

        self._set_sources = False

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

        assert not self._set_sources

        self._set_sources = True

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

#pylint: enable=too-few-public-methods
