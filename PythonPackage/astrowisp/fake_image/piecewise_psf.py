"""Declares the PiecewisePSF base class for piecewise PSF functions."""

from bisect import bisect

from astrowisp.psf_base import PSFBase

def get_piece_index(boundaries, coordinate):
    """
    Return the index of the piece along one axis containing the given coord.

    Args:
        boundaries:    The offsets relative to the PSF center where different
            PSF pieces meet along the direction in which we are trying to locate
            the piece index.
        coordinate:    The coordinate we are trying to find the piece index of.

    Returns:
        ind:    The index along the selected coordinate of the piece
            containing x.
    """

    if coordinate == boundaries[-1]:
        return len(boundaries) - 2
    return bisect(boundaries, coordinate) - 1

class PiecewisePSF(PSFBase):
    """Base clas for PSFs defined on a grid of pieces."""

    def __init__(self, boundaries, pieces):
        """
        Define a PSF with the given boundaries.

        Args:
            boundaries:    Dictionary with keys `x` and `y` giving the offsets
                relative to the center of the horizontal piece boundaries. The
                PSF is zero left of the first or right of the last x boundary
                as wall as below the first and above the last y boundary.

            pieces:    The pieces making up the PSF should be a class inherited
                from PSFPiece.

        Returns: None
        """

        self._boundaries = boundaries
        self._pieces = pieces
        super().__init__()

    def get_left_range(self):
        """Return how far the PSF extends to the left of center."""

        return -self._boundaries['x'][0]

    def get_right_range(self):
        """Return how far the PSF extends to the right of center."""

        return self._boundaries['x'][-1]

    def get_down_range(self):
        """Return how far the PSF extends downward of center."""

        return -self._boundaries['y'][0]

    def get_up_range(self):
        """Return how far the PSF extends upward of center."""

        return self._boundaries['y'][-1]

    #(x, y) is a reasonable way to specify the coordinates of an offset vector.
    #pylint: disable=invalid-name
    def __call__(self, x, y):
        """See the documentation of PSFBase.__call__."""

        if(
                x < self._boundaries['x'][0]
                or
                x > self._boundaries['x'][-1]
                or
                y < self._boundaries['y'][0]
                or
                y > self._boundaries['y'][-1]
        ):
            return 0.0

        return self._pieces[
            get_piece_index(self._boundaries['y'], y)
        ][
            get_piece_index(self._boundaries['x'], x)
        ](x, y)
    #pylint: enable=invalid-name

    def integrate(self, left, bottom, width, height):
        """See documentation of PSFBase.integrate."""

        result = 0.0

        if width < 0:
            sign = -1
            left = left + width
            width = abs(width)
        else:
            sign = 1
        if height < 0:
            sign *= -1
            bottom = bottom + height
            height = abs(height)

        for y_piece in range(
                max(0, get_piece_index(self._boundaries['y'], bottom)),
                min(get_piece_index(self._boundaries['y'],
                                    bottom + height) + 1,
                    len(self._boundaries['y']) - 1)
        ):

            y_min = max(self._boundaries['y'][y_piece], bottom)
            y_max = min(height + bottom, self._boundaries['y'][y_piece + 1])

            for x_piece in range(
                    max(0, get_piece_index(self._boundaries['x'], left)),
                    min(get_piece_index(self._boundaries['x'],
                                        left + width) + 1,
                        len(self._boundaries['x']) - 1)
            ):
                piece = self._pieces[y_piece][x_piece]
                x_min = max(self._boundaries['x'][x_piece], left)
                x_max = min(width + left, self._boundaries['x'][x_piece + 1])

                result += piece.integrate(x_min,
                                          y_min,
                                          x_max - x_min,
                                          y_max - y_min)
        return sign * result
