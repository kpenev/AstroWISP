from fake_image.PSFPiece import PSFPiece

class BipolynomialPSFPiece(PSFPiece) :
    """Class for PSF pieces over which the PSF is a bi-cubic function."""

    def __init__(self, coefficients) :
        """
        Initialize a bi-cubic PSF piece with a set of coefficients.

        Args:
            - coefficients:
                The bi-cubic polynomial coefficients. Should be a two-index
                structure with the first index beige y power of the term and
                the second index the x-power.

        Returns: None
        """

        self.__coefficients = coefficients

    def __call__(self, x, y) :
        """Evaluate the bi-cubic polynomial at the given position."""

        result = 0.0
        y_factor = 1.0
        for y_pow in range(len(self.__coefficients)) :
            x_factor = 1.0
            for x_pow in range(len(self.__coefficients[y_pow])) :
                result += (self.__coefficients[y_pow][x_pow]
                           *
                           x_factor
                           *
                           y_factor)
                x_factor *= x
            y_factor *= y

        return result

    def integrate(self, left, bottom, width, height) :
        """
        Evaluate the integral of the bi-cubic polynomial over a rectangle.

        Args:
            - left:
                The x-coordinate of the left boundary of the rectangle to
                integrate over.

            - bottom:
                The y-coordinate of the bottom boundary of the rectangle to
                integrate over.

            - width:
                The x-size of the rectangle to integrate over.

            -height:
                The y-size of the rectangle to integrate over.

        Returns:
            The integral of the bi-cubic polynomial function defining this
            piece of the PSF over the specified rectangle, with no
            consideration of whether the rectangle fits within the piece.
        """

        result = 0.0
        top = bottom + height
        right = left + width
        bottom_factor = bottom
        top_factor = top
        for y_pow in range(len(self.__coefficients)) :
            left_factor = left
            right_factor = right
            for x_pow in range(len(self.__coefficients[y_pow])) :
                result += (
                    self.__coefficients[y_pow][x_pow]
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
