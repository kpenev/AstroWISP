"""Define the BipolynomialPSFPiece class."""

from astrowisp.fake_image.psf_piece import PSFPiece

#This is only a base class, so eventual objects will have more public methods
#pylint: disable=too-few-public-methods
class BipolynomialPSFPiece(PSFPiece):
    """Class for PSF pieces over which the PSF is a bi-polynomial function."""

    def __init__(self, coefficients):
        """
        Initialize a bi-cubic PSF piece with a set of coefficients.

        Args:
            coefficients:    The bi-cubic polynomial coefficients. Should be a
                two-index structure with the first index beige y power of the
                term and the second index the x-power.

        Returns:
            None
        """

        self.__coefficients = coefficients
        super().__init__()

    #(x, y) is a perfectly valid way to specify a point.
    #pylint: disable=invalid-name
    def __call__(self, x, y):
        """Evaluate the bi-cubic polynomial at the given position."""

        result = 0.0
        y_factor = 1.0
        for ypow_coef_slice in self.__coefficients:
            x_factor = 1.0
            for coef in ypow_coef_slice:
                result += (coef
                           *
                           x_factor
                           *
                           y_factor)
                x_factor *= x
            y_factor *= y

        return result
    #pylint: enable=invalid-name

    #It's either that or not using enumarete
    #pylint: disable=too-many-locals
    def integrate(self, left, bottom, width, height):
        """See documentation of PSFPiece.integrate."""

        result = 0.0
        top = bottom + height
        right = left + width
        bottom_factor = bottom
        top_factor = top
        for y_pow, ypow_coef_slice in enumerate(self.__coefficients):
            left_factor = left
            right_factor = right
            for x_pow, coef in enumerate(ypow_coef_slice):
                result += (
                    coef
                    *
                    (right_factor - left_factor) / (x_pow + 1)
                    *
                    (top_factor - bottom_factor) / (y_pow + 1)
                )
                left_factor *= left
                right_factor *= right
            bottom_factor *= bottom
            top_factor *= top

        return result
    #pylint: enable=too-many-locals
#pylint: enable=too-few-public-methods
