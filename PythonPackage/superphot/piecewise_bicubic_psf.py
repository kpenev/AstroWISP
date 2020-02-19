"""A wrapper class for working with PSFs/PRFs from the C/C++ library."""

from ctypes import c_double
import numpy

from superphot.psf_base import PSFBase
from superphot._initialize_library import superphot_library

class PiecewiseBicubicPSF(PSFBase):
    """Implement the PSFBase methods for libary PSFs."""

    def __init__(self, library_psf):
        """Wrap the given library PSF in a convenient python interface."""

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

        print('Evaluating PRF at x: ' + repr(x))
        print('Evaluating PRF at y: ' + repr(y))
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

        superphot_library.evaluate_piecewise_bicubic_psf(
            self._library_psf,
            x,
            y,
            x.size,
            result
        )

        return result
    #pylint: enable=invalid-name

    def __del__(self):
        """Delete the underlying library PSF."""

        superphot_library.destroy_piecewise_bicubic_psf(self._library_psf)

    def get_left_range(self):
        raise NotImplementedError

    def get_right_range(self):
        raise NotImplementedError

    def get_down_range(self):
        raise NotImplementedError

    def get_up_range(self):
        raise NotImplementedError

    def integrate(self, left, bottom, width, height):
        raise NotImplementedError
