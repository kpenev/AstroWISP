"""Defines a base class for all PSF models: PSFBase."""

from abc import ABC, abstractmethod

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
    def integrate(self, left, bottom, width, height):
        """
        Return the integral of the PSF over a rectangle.

        Args:
            left:    The horizontal offset relative to the center of the left
                boundary of the rectangle to integrate over.

            bottom:    The vertical offset relative to the center of the bottom
                boundary of the rectangle to integrate over.

            width:    The horizontal size of the rectangle to integrate over.

            height:    The vertical size of the rectangle to integrate over.

        Returns:
            integral:    The integral of the PSF over left < x < left + width
                bottom < y < bottom + height
        """
