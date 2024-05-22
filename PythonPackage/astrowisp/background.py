"""Interface to the AstroWISP background library."""

from ctypes import c_double, c_uint
import numpy

from astrowisp._initialize_library import get_astrowisp_library

#The __init__, __del__ and __call__ methods justify making this a class.
#pylint: disable=too-few-public-methods
#pylint: disable=too-many-instance-attributes
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

        self._astrowisp_library = get_astrowisp_library()
        self.image = image
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.error_confidence = error_confidence
        self._library_image = self._astrowisp_library.create_core_image(
            self.image.shape[1],
            self.image.shape[0],
            self.image,
            None,
            None,
            True
        )
        self.library_extractor = (
            self._astrowisp_library.create_background_extractor(
                inner_radius,
                outer_radius,
                inner_radius,
                self._library_image,
                error_confidence
            )
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

        self._astrowisp_library.add_source_list_to_background_extractor(
            self.library_extractor,
            source_x,
            source_y,
            source_x.size
        )

        bg_value = numpy.empty(source_x.size, dtype=c_double)
        bg_error = numpy.empty(source_x.size, dtype=c_double)
        bg_numpix = numpy.empty(source_x.size, dtype=c_uint)
        self._astrowisp_library.get_all_backgrounds(
            self.library_extractor,
            bg_value,
            bg_error,
            bg_numpix
        )
        return bg_value, bg_error, bg_numpix

    def __del__(self):
        r"""Destroy the image and extractor created in :meth:`__init__`\ ."""

        self._astrowisp_library.destroy_core_image(self._library_image)
        self._astrowisp_library.destroy_background_extractor(
            self.library_extractor
        )
#pylint: enable=too-few-public-methods
#pylint: enable=too-many-instance-attributes
