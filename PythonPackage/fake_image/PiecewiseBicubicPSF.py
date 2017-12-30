from fake_image.PiecewisePSF import PiecewisePSF
from fake_image.BipolynomialPSFPiece import BipolynomialPSFPiece
import scipy.linalg

class PiecewiseBicubicPSF(PiecewisePSF) :
    """A PSF class where the PSF over each piece is a bi-cubic function."""

    def _create_piece(self,
                      x_boundaries,
                      y_boundaries,
                      values,
                      d_dx,
                      d_dy,
                      d2_dxdy) :
        """
        Return a BicubicPSFPiece satisfying the given constraints.

        Args:
            - x_boundaries:
                A size=2 structure with the first/second entry being the x
                coordinate of the left/right boundary of the cell.

            - y_boundaries:
                A size=2 structure with the first/second entry being the y
                coordinate of the bottom/top boundary of the cell.

            - values:
                A 2x2 structure giving the values of the piece bi-cubic
                polynomial af the corners of the piece edge.

            - d_dx:
                A 2x2 structure giving the x derivatives of the piece
                bi-cubic polynomial af the corners of the piece edge.

            - d_dy:
                A 2x2 structure giving the y derivatives of the piece
                bi-cubic polynomial af the corners of the piece edge.

            - d2_dxdy:
                A 2x2 structure giving the x,y cross-derivatives of the piece
                bi-cubic polynomial af the corners of the piece edge.

        Returns:
            A BicubicPSFPiece instance with value, x, y and xy derivatives of
            values[i][j], d_dx[i][j], d_dy[i][j], d2_dxdy[i][j] at
            (x_boundaries[i], y_boundaries[j]).
        """

        matrix = scipy.empty((16,16))
        row_offset = 0
        for vert_index in range(2) :
            y = y_boundaries[vert_index]
            for horz_index in range(2) :
                x = x_boundaries[horz_index]
                y_term = 1.0
                column = 0
                for y_pow in range(4) :
                    x_term = 1.0
                    for x_pow in range(4) :
                        if x == 0 :
                            dx_term = (1.0 if x_pow == 1 else 0.0)
                        else :
                            dx_term = x_pow * x_term / x

                        if y == 0 :
                            dy_term = (1.0 if y_pow == 1 else 0.0)
                        else :
                            dy_term = y_pow * y_term / y

                        matrix[row_offset, column] = x_term * y_term
                        matrix[row_offset + 4, column] = dx_term * y_term
                        matrix[row_offset + 8, column] = x_term * dy_term
                        matrix[row_offset + 12, column] = dx_term * dy_term

                        column += 1
                        x_term *= x
                    y_term *= y
                row_offset += 1

        rhs = scipy.empty(16)
        rhs[ 0 : 4 ] = values.flatten()
        rhs[ 4 : 8 ] = d_dx.flatten()
        rhs[ 8 : 12 ] = d_dy.flatten()
        rhs[ 12 : 16 ] = d2_dxdy.flatten()
        coefficients = scipy.linalg.solve(matrix, rhs)
        return BipolynomialPSFPiece(coefficients.reshape((4,4)))

    def __init__(self,
                 x_boundaries,
                 y_boundaries,
                 values,
                 d_dx,
                 d_dy,
                 d2_dxdy) :
        """
        Initialize a PiecewiseBicubicPSF with the given shape.

        Args:
            - x_boundaries:
                See PiecewisePSF.__init__

            - y_boundaries:
                See PiecewisePSF.__init__

            - values:
                The value of the PSF at each of the grid intersections.
                Should be a 2-index object with the first index iterating
                over y boundaries and the second index over x_boundaries.

            - d_dx:
                The derivative w.r.t. x of the PSF at each of the grid
                intersections. Same format as values.

            - d_dy:
                The derivative w.r.t. y of the PSF at each of the grid
                intersections. Same format as values.

            - d2_dxdy:
                The mixed derivative w.r.t. x and y of the PSF at each of the
                grid intersections. Same format as values.

        Returns: None
        """
        
        pieces = []
        for cell_y_index in range(len(y_boundaries) - 1) :
            pieces.append([])
            for cell_x_index in range(len(x_boundaries) - 1) :
                pieces[-1].append(
                    self._create_piece(
                        x_boundaries[cell_x_index : cell_x_index + 2],
                        y_boundaries[cell_y_index : cell_y_index + 2],
                        values[cell_y_index : cell_y_index + 2,
                               cell_x_index : cell_x_index + 2],
                        d_dx[cell_y_index : cell_y_index + 2,
                             cell_x_index : cell_x_index + 2],
                        d_dy[cell_y_index : cell_y_index + 2,
                             cell_x_index : cell_x_index + 2],
                        d2_dxdy[cell_y_index : cell_y_index + 2,
                                cell_x_index : cell_x_index + 2]
                    )
                )
        super().__init__(x_boundaries, y_boundaries, pieces)
