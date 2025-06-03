"""Define the :class:`SubPixPhot`, which performs aperture photometry."""

from ctypes import c_double, c_uint

import numpy

from astrowisp._initialize_library import get_astrowisp_library


# This only makes sense as a class.
# pylint: disable=too-few-public-methods
class SubPixPhot:
    """
    Use sub-pixel aware aperture photomotre to measure fluxes of sources.

    Attributes:
        _library_configuration (_c_subpixphot_configuration):    The library
            configuration object to use the next aperture photometry
            measurement.

        _result_tree (_c_h5io_data_tree_p):    The IOTree instance
            containing PSF fitting and aperture photometry information from the
            last aperture photomtry measurement, or None if no photometry has
            been done yet.

        _library_subpix_map (_c_core_sub_pixel_map_p):    The library object
            representing the currently set sub-pixel map.

        configuration (dict):    The congfiguration for how to perform the
            next aperture photometry. The following keys are used (others are
            ignored by this class).

            subpixmap (2D numpy array):
                The sub-pixel map to assume.

            apertures ([float, ...]):
                A list of the apertures to use.

            gain (float):
                The gain to assume for the input image.

            magnitude_1adu (float):
                The magnitude that corresponds to a flux of 1ADU.

            const_error (float):
                A constant that gets added to all error estimates.

        image (dict):    The last image for which aperture photometry was
            extracted. Contains the following entries:

                values:
                    The calibrated pixel values.

                errors:
                    Error estimates for the entries in values.

                mask:
                    Bitmask flags indicating any known problems with image
                    pixels.
    """

    _default_configuration = {
        "subpixmap": numpy.ones((1, 1), dtype=c_double),
        "apertures": numpy.arange(1.0, 5.5),
        "gain": 1.0,
        "magnitude_1adu": 10.0,
        "const_error": 0.0,
    }

    @staticmethod
    def _format_config(param_value):
        """
        Format config param for passing to AstroWISP aperture photometry func.

        Args:
            param_value (pair of values):    One of the recognized keys from
                :attr:`configuration` and the value it should be set to.

        Returns:
            (bytes, bytes):
                * The name of the option corresponding to the configuration
                  parameter being set.
                * an ascii string representing the value
                  in the format expected by the configuration file parser of the
                  subpixphot tool.
        """

        if param_value[0] == "subpixmap":
            return ()

        if param_value[0] == "apertures":
            return (
                b"ap.aperture",
                b",".join([repr(ap).encode("ascii") for ap in param_value[1]]),
            )

        if param_value[0] == "const_error":
            param_name = b"ap.const-error"
        else:
            param_name = param_value[0].replace("_", "-").encode("ascii")

        return (param_name, repr(param_value[1]).encode("ascii"))

    def __init__(self, **configuration):
        r"""
        Prepare an object for measuring fluxes using aperture photometry.

        Args:
            **configuration:    See :attr:`configuration`\ .

        Returns:
            None
        """

        self._astrowisp_library = get_astrowisp_library()
        self.image = None
        self._library_subpix_map = None

        self._library_configuration = (
            self._astrowisp_library.create_subpixphot_configuration()
        )

        self.configuration = dict(self._default_configuration)
        self.configuration.update(**configuration)
        self.configure(**self.configuration)
        self._result_tree = None

    def configure(self, **configuration):
        r"""
        Modify the currently defined configuration.

        Args:
            **configuration:    See :attr:`configuration`\ .

        Returns:
            None
        """

        for k, value in configuration.items():
            if k not in self.configuration:
                raise KeyError(
                    "Unrecognized configuration parameter: " + repr(k)
                )
            if k == "subpixmap":
                if self._library_subpix_map is not None:
                    self._astrowisp_library.destroy_core_subpixel_map(
                        self._library_subpix_map
                    )
                self._library_subpix_map = (
                    self._astrowisp_library.create_core_subpixel_map(
                        *value.shape, value
                    )
                )

        self.configuration.update(configuration)

        config_arguments = sum(
            map(self._format_config, self.configuration.items()), ()
        ) + (b"",)

        self._astrowisp_library.update_subpixphot_configuration(
            self._library_configuration, *config_arguments
        )

    # No clean way to reduce the number of argumets.
    # pylint: disable=too-many-arguments
    def __call__(self, image, fitpsf_io_tree, image_index=0):
        r"""
        Measure the fluxes of all sources in an image using aperture photometry.

        Args:
            image(2-D numpy.array, 2-D numpy.array, 2-D numpy.array): entries:

                1. The calibrated values of the pixels in the image.

                2. Error estimates of the pixel values

                3. Bitmask flagging any known issues with image pixels (e.g.
                   saturation, hot pixels, etc.).

            fitsf_io_tree (IOTree):    The result tree returned by
                fit_star_shape.fit(). On output, this variable also contains the
                newly derived aperture photometry measurements.

            image_index (int):    The image index within the result tree
                corresponding to the input image. This is the index of the image
                within the list of images passed to :class:`fit_star_shape`\ .

        Returns:
            None. However, the ``fitpsf_io_tree`` argument is updated to
            include the newly measured fluxes.
        """

        assert image[1].shape == image[0].shape
        assert image[2].shape == image[0].shape

        self.image = image

        library_image = self._astrowisp_library.create_core_image(
            image[0].shape[1], image[0].shape[0], *image, True
        )

        self._astrowisp_library.subpixphot(
            library_image,
            self._library_subpix_map,
            self._library_configuration,
            fitpsf_io_tree.library_tree,
            c_uint(image_index),
        )

        self._astrowisp_library.destroy_core_image(library_image)

    # pylint: enable=too-many-arguments

    def __del__(self):
        r"""Destroy any library objects created by this object."""

        self._astrowisp_library.destroy_subpixphot_configuration(
            self._library_configuration
        )
        if self._library_subpix_map is not None:
            self._astrowisp_library.destroy_core_subpixel_map(self._library_subpix_map)


# pylint: enable=too-few-public-methods
