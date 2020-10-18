"""Defines a base class for all PSF models: PSFBase."""

from abc import ABC, abstractmethod

import numpy

class PSFBase(ABC):
    """The base class for all supported PSFs."""

    @abstractmethod
    #(x, y) is a reasonable way to specify the components of an offset vector.
    #pylint: disable=invalid-name
    def __call__(self, x, y):
        """
        Return the value of the PSF at the given point.

        Args:
            x:    The horizontal offset from center of the point at which to
                return the PSF value.

            y:    The vertical offset from center of the point at which to
                return the PSF value.

        Returns:
            The value of the PSF at (x, y) relative to the source center.
        """
    #pylint: enable=invalid-name

    @abstractmethod
    def get_left_range(self):
        """Return how far the PSF extends to the left of center."""

    @abstractmethod
    def get_right_range(self):
        """Return how far the PSF extends to the right of center."""

    @abstractmethod
    def get_down_range(self):
        """Return how far the PSF extends downward of center."""

    @abstractmethod
    def get_up_range(self):
        """Return how far the PSF extends upward of center."""

    @abstractmethod
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

    def predict_pixel(self,
                      center_x,
                      center_y,
                      subpix_map=numpy.array([[1.0]])):
        """
        Predict the value(s) of pixel(s) assuming a pixel sensitivity map.

        Args:
            center_x(float or array):    The offset(s) of the center(s) of the
                pixels from the source position in the x direction.

            center_y(float or array):    The offset(s) of the center(s) of the
                pixels from the source position in the y direction.

            subpix_map(2-D array):    The pixel sensitivity map to assume.

        Returns:
            float or array:
                The integral(s) of product of the PSF and the sub-pixel
                sensitivity map over pixel(s) located as specified relative to
                the source center.
        """

        if isinstance(center_x, numpy.ndarray):
            if not isinstance(center_y, numpy.ndarray):
                center_y = numpy.full(center_x.shape, center_y)
        else:
            if isinstance(center_y, numpy.ndarray):
                center_x = numpy.full(center_y.shape, center_x)
            else:
                center_x = numpy.array([float(center_x)])
                center_y = numpy.array([float(center_y)])


        width = 1.0 / subpix_map.shape[1]
        height = 1.0 / subpix_map.shape[0]
        subpix_center_x, subpix_center_y = numpy.meshgrid(
            numpy.arange(-0.5, 0.5, width) + width / 2.0,
            numpy.arange(-0.5, 0.5, height) + height / 2.0
        )
        subpix_center_x = center_x.ravel()[:, None] + subpix_center_x.ravel()
        subpix_center_y = center_y.ravel()[:, None] + subpix_center_y.ravel()
        return self.integrate(
            center_x=subpix_center_x,
            center_y=subpix_center_y,
            width=numpy.full(subpix_center_x.shape, width),
            height=numpy.full(subpix_center_x.shape, height)
        ).dot(
            subpix_map.ravel()
        )
