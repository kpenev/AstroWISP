"""A wrapper class for working with PSFs/PRFs from the C/C++ library."""

from ctypes import c_double
import numpy

from astrowisp.psf_base import PSFBase
from astrowisp._initialize_library import get_astrowisp_library

class PiecewiseBicubicPSF(PSFBase):
    """Implement the PSFBase methods for libary PSFs."""

    def __init__(self, library_psf):
        """Wrap the given library PSF in a convenient python interface."""

        self._astrowisp_library = get_astrowisp_library()
        self._library_psf = library_psf

    #(x, y) is a reasonable way to specify the components of an offset vector.
    #pylint: disable=invalid-name
    def __call__(self, x, y):
        """
        Return the value(s) of the PSF at the given point(s).

        Args:
            x(float or numpy array):    The horizontal offset(s) from center of
                the point at which to return the PSF value.

            y(float or numpy array):    The vertical offset(s) from center of
                the point at which to return the PSF value.

        Returns:
            numpy array:
                The value(s) of the PSF at (x, y) relative to the source center.
        """

        print('Evaluating PRF at x with shape '
              +
              repr(numpy.atleast_1d(x).shape))
        print('Evaluating PRF at y with shape '
              +
              repr(numpy.atleast_1d(y).shape))
        if isinstance(x, numpy.ndarray):
            if not isinstance(y, numpy.ndarray):
                y = numpy.full(x.shape, y)
        else:
            if isinstance(y, numpy.ndarray):
                x = numpy.full(y.shape, x)
            else:
                x = numpy.array([float(x)])
                y = numpy.array([float(y)])

        assert x.size == y.size

        result = numpy.empty(x.size, dtype=c_double)

        self._astrowisp_library.evaluate_piecewise_bicubic_psf(
            self._library_psf,
            x.flatten(),
            y.flatten(),
            x.size,
            result
        )

        return result.reshape(x.shape)
    #pylint: enable=invalid-name

    def integrate(self,
                  *,
                  center_x,
                  center_y,
                  width,
                  height,
                  circle_radius=None):
        """
        Return integrals of the PSF over circle-rectangle overlaps.

        Args:
            center_x(float or array):    The x coordinate(s) of the center(s) of
                the rectangle(s) to integrate over.

            center_y(float or array):    The y coordinate(s) of the center(s) of
                the rectangle(s) to integrate over.

            width(float or array):    The width(s) of the rectangle(s).

            height(float or array):    The height(s) of the rectangle(s).

            circle_radius(float or array):    The rad(ii/us) of the circle(s).
                For zero entries or None, the integral is over the full
                rectangle(s).

        Returns:
            float or array:
                The integral of the PSF over the specified area(s).
        """

        if not isinstance(center_x, numpy.ndarray):
            center_x = numpy.array([center_x])
            center_y = numpy.array([center_y])
            width = numpy.array([width])
            height = numpy.array([height])
            circle_radius = numpy.array(
                [0.0 if circle_radius is None else circle_radius]
            )
        else:
            assert center_x.shape == center_y.shape
            assert center_x.shape == width.shape
            assert center_x.shape == height.shape

        if circle_radius is None:
            circle_radius = numpy.full(center_x.shape, 0.0)
        else:
            assert center_x.shape == circle_radius.shape

        result = numpy.empty(center_x.shape, dtype=float, order='C')

        self._astrowisp_library.integrate_piecewise_bicubic_psf(
            self._library_psf,
            center_x.astype(c_double, order='C', copy=False).ravel(),
            center_y.astype(c_double, order='C', copy=False).ravel(),
            width.astype(c_double, order='C', copy=False).ravel(),
            height.astype(c_double, order='C', copy=False).ravel(),
            circle_radius.astype(c_double, order='C', copy=False).ravel(),
            result.size,
            result.astype(c_double, order='C', copy=False).ravel()
        )

        if not isinstance(center_x, numpy.ndarray):
            result = float(result)

        return result

    def __del__(self):
        """Delete the underlying library PSF."""

        self._astrowisp_library.destroy_piecewise_bicubic_psf(self._library_psf)

    def get_left_range(self):
        raise NotImplementedError

    def get_right_range(self):
        raise NotImplementedError

    def get_down_range(self):
        raise NotImplementedError

    def get_up_range(self):
        raise NotImplementedError
