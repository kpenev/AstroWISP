"""Defines the PiecewiseBicubicPSF class."""

import scipy.linalg
import numpy

from astrowisp.fake_image.piecewise_psf import PiecewisePSF
from astrowisp.fake_image.bipolynomial_psf_piece import BipolynomialPSFPiece

class PiecewiseBicubicPSF(PiecewisePSF):
    """
    A piecewise PSF class where the PSF over each piece is a bi-cubic function.
    """

    #Could not think of a reasonable way to split this function or decrease the
    #number of local variables without losing readability.
    #pylint: disable=too-many-locals
    @staticmethod
    def _create_piece(boundaries, psf_parameters):
        """
        Return a BicubicPSFPiece satisfying the given constraints.

        Args:
            The same as the arguments of __init__

        Returns:
            psf_piece:
                A BicubicPSFPiece instance with parameters as specified by
                `psf_parameters`.
        """

        matrix = numpy.empty((16, 16))
        row_offset = 0
        for vert_index in range(2):
            #y is a reasonable name for position vector y-component.
            #pylint: disable=invalid-name
            y = boundaries['y'][vert_index]
            #pylint: enable=invalid-name
            for horz_index in range(2):
                #y is a reasonable name for position vector y-component.
                #pylint: disable=invalid-name
                x = boundaries['x'][horz_index]
                #pylint: enable=invalid-name
                y_term = 1.0
                column = 0
                for y_pow in range(4):
                    x_term = 1.0
                    for x_pow in range(4):
                        if x == 0:
                            dx_term = (1.0 if x_pow == 1 else 0.0)
                        else:
                            dx_term = x_pow * x_term / x

                        if y == 0:
                            dy_term = (1.0 if y_pow == 1 else 0.0)
                        else:
                            dy_term = y_pow * y_term / y

                        matrix[row_offset, column] = x_term * y_term
                        matrix[row_offset + 4, column] = dx_term * y_term
                        matrix[row_offset + 8, column] = x_term * dy_term
                        matrix[row_offset + 12, column] = dx_term * dy_term

                        column += 1
                        x_term *= x
                    y_term *= y
                row_offset += 1

        rhs = numpy.empty(16)
        rhs[0 : 4] = psf_parameters['values'].flatten()
        rhs[4 : 8] = psf_parameters['d_dx'].flatten()
        rhs[8 : 12] = psf_parameters['d_dy'].flatten()
        rhs[12 : 16] = psf_parameters['d2_dxdy'].flatten()
        coefficients = scipy.linalg.solve(matrix, rhs)
        return BipolynomialPSFPiece(coefficients.reshape((4, 4)))
    #pylint: enable=too-many-locals

    def __init__(self,
                 boundaries,
                 psf_parameters):
        """
        Initialize a PiecewiseBicubicPSF with the given shape.

        Args:
            boundaries:    Dictionary (keys ``x`` and ``y``) listing the cell
                horizontal/vertical boundaries.

            psf_parameters:    A dictionary of 2x2 structures with keys:

                * values:    The values of the piece bi-cubic polynomial af the
                      intersections of the horizontal & vertical ``boundaries``.

                * d_dx:    The x derivatives of the piece bi-cubic polynomial af
                      the intersections of the horizontal & vertical
                      ``boundaries``.

                * d_dy:    The y derivatives of the piece bi-cubic polynomial af
                      the intersections of the horizontal & vertical
                      ``boundaries``.

                * d2_dxdy:    The x,y cross-derivatives of the piece bi-cubic
                      polynomial af the intersections of the horizontal &
                      vertical ``boundaries``.

        Returns:
            None
        """

        pieces = []
        for cell_y_index in range(len(boundaries['y']) - 1):
            pieces.append([])
            for cell_x_index in range(len(boundaries['x']) - 1):
                pieces[-1].append(
                    self._create_piece(
                        boundaries={
                            'x': boundaries['x'][cell_x_index
                                                 :
                                                 cell_x_index + 2],
                            'y': boundaries['y'][cell_y_index
                                                 :
                                                 cell_y_index + 2]
                        },
                        psf_parameters={
                            'values': psf_parameters['values'][
                                cell_y_index : cell_y_index + 2,
                                cell_x_index : cell_x_index + 2
                            ],
                            'd_dx': psf_parameters['d_dx'][
                                cell_y_index : cell_y_index + 2,
                                cell_x_index : cell_x_index + 2
                            ],
                            'd_dy': psf_parameters['d_dy'][
                                cell_y_index : cell_y_index + 2,
                                cell_x_index : cell_x_index + 2
                            ],
                            'd2_dxdy': psf_parameters['d2_dxdy'][
                                cell_y_index : cell_y_index + 2,
                                cell_x_index : cell_x_index + 2
                            ]
                        }
                    )
                )
        super().__init__(boundaries, pieces)
