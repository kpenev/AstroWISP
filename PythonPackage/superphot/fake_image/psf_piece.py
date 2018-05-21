"""Define the base class for pieces of piecewise PSFs."""

from abc import ABC, abstractmethod

#This is only a base class, so eventual objects will have more public methods
#pylint: disable=too-few-public-methods
class PSFPiece(ABC):
    """Declare a minimum interface for pieces of PiecewisePSF."""

    #(x,y) is a reasonable way to specify the components of an offset vector
    #pylint: disable=invalid-name
    @abstractmethod
    def __call__(self, x, y):
        """Evaluate the piece at the given position."""

    @abstractmethod
    def integrate(self, left, bottom, width, height):
        """
        Evaluate the integral of the cell function over a rectangle.

        Args:
            left:    The x-coordinate of the left boundary of the rectangle to
                integrate over.

            bottom:    The y-coordinate of the bottom boundary of the
                rectangle to integrate over.

            width:    The x-size of the rectangle to integrate over.

            height:    The y-size of the rectangle to integrate over.

        Returns:
            float:
                The integral of the bi-cubic polynomial function defining this
                piece of the PSF over the specified rectangle, with no
                consideration of whether the rectangle fits within the piece.
        """
    #pylint: enable=invalid-name

#pylint: enable=too-few-public-methods
