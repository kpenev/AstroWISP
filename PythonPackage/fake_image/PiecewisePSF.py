from fake_image.PSFBase import PSFBase
from math import floor, ceil
from bisect import bisect

def get_piece_index(boundaries, x) :
    """
    Return the index of the piece along one axis containing the given coord.
    
    Args:
        - boundaries:
            The offsets relative to the PSF center where different PSF pieces
            meet along the direction in which we are trying to locate the
            piece index.
        - x:
            The coordinate we are trying to find the piece index of.

    Returns:
        The index along the selected coordinate of the piece containing x.
    """

    if x == boundaries[-1] : return len(boundaries) - 2
    return bisect(boundaries, x) - 1

class PiecewisePSF(PSFBase) :
    """
    Base clas for PSFs defined on a grid of pieces.
    """

    def __init__(self, x_boundaries, y_boundaries, pieces) :
        """
        Define a PSF with the given boundaries.

            - x_boundaries:
                The offsets relative to the center of the horizontal piece
                boundaries. The PSF is zero left of the first or right of the
                last boundary given.

            - y_boundaries:
                The offsets relative to the center of the vertical piece
                boundaries. The PSF is zero downward of the first or upward
                of the last boundary given.

        Returns: None
        """

        self.__x_boundaries = x_boundaries
        self.__y_boundaries = y_boundaries
        self._pieces = pieces

    def get_left_range(self) :
        """Return how far the PSF extends to the left of center."""

        return -self.__x_boundaries[0]

    def get_right_range(self) :
        """Return how far the PSF extends to the right of center."""

        return self.__x_boundaries[-1]

    def get_down_range(self) :
        """Return how far the PSF extends downward of center."""

        return -self.__y_boundaries[0]

    def get_up_range(self) :
        """Return how far the PSF extends upward of center."""

        return self.__y_boundaries[-1]

    def __call__(self, x, y) :
        """
        Return the value of the PSF at the given point.

        Args:
            - x:
                The horizontal offset from center of the point at which to
                return the PSF value.

            - y:
                The vertical offset from center of the point at which to
                return the PSF value.

        Returns:
            The value of the PSF at (x, y) relative to the source center.
        """

        if(
                x < self.__x_boundaries[0]
                or
                x > self.__x_boundaries[-1]
                or
                y < self.__y_boundaries[0]
                or
                y > self.__y_boundaries[-1]
        ) :
            return 0.0
        else :
            return self._pieces[
                get_piece_index(self.__y_boundaries, y)
            ][
                get_piece_index(self.__x_boundaries, x)
            ](x, y)

    def integrate(self, left, bottom, width, height) :
        """
        Return the integral of the PSF over a rectangle.

        Args:
            - left:
                The horizontal offset relative to the center of the left
                boundary of the rectangle to integrate over.

            - bottom:
                The vertical offset relative to the center of the bottom
                boundary of the rectangle to integrate over.

            - width:
                The horizontal size of the rectangle to integrate over.

            - height:
                The vertical size of the rectangle to integrate over.

        Returns:
            The integral of the PSF over
                left < x < left + width
                bottom < y < bottom + height
        """

        result = 0.0

        if width < 0 :
            sign = -1
            left = left + width
            width = abs(width)
        else :
            sign = 1
        if height < 0 :
            sign *= -1
            bottom = bottom + height
            height = abs(height)

        for y_piece in range(
                max(0, get_piece_index(self.__y_boundaries, bottom)),
                min(get_piece_index(self.__y_boundaries,
                                    bottom + height) + 1,
                    len(self.__y_boundaries) - 1)
        ) :

            y_min = max(self.__y_boundaries[y_piece], bottom)
            y_max = min(height + bottom, self.__y_boundaries[y_piece + 1])

            for x_piece in range(
                    max(0, get_piece_index(self.__x_boundaries, left)),
                    min(get_piece_index(self.__x_boundaries,
                                        left + width) + 1,
                        len(self.__x_boundaries) - 1)
            ) :
                piece = self._pieces[y_piece][x_piece]
                x_min = max(self.__x_boundaries[x_piece], left)
                x_max = min(width + left, self.__x_boundaries[x_piece + 1])

                result += piece.integrate(x_min,
                                          y_min,
                                          x_max - x_min,
                                          y_max - y_min)
        return sign * result
