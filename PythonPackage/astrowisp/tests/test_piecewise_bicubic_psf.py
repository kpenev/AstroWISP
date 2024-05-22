#!/usr/bin/env python3

"""Define unittest test case for the PiecewiseBicubicPSF class."""

import unittest
import numpy

from astrowisp.tests.utilities import FloatTestCase
from astrowisp.fake_image.piecewise_bicubic_psf import PiecewiseBicubicPSF


def point_in_grid(point, grid):
    """Return True iff the point is within grid."""

    return (
        point['x'] >= grid['x'][0]
        and
        point['x'] <= grid['x'][-1]
        and
        point['y'] >= grid['y'][0]
        and
        point['y'] <= grid['y'][-1]
    )

class TestPiecewiseBicubicPSF(FloatTestCase):
    """Make sure the PiecewiseBicubicPSF class functions as expected."""

    def setUp(self):
        """Define a set of picewise PSF grids to use during tests."""

        self._grids = [
            {'x': numpy.array([0.0, 1.0]),
             'y': numpy.array([0.0, 1.0])},
            {'x': numpy.array([0.0, 0.5]),
             'y': numpy.array([0.0, numpy.pi / 3.0])},
            {'x': numpy.array([0.0, 0.5, 1.0]),
             'y': numpy.array([0.0, 1.0])},
        ]
        self._test_points = [
            [
                {'x': x, 'y': y}
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2,
                                        numpy.pi, 10)
            ],
            [
                {'x': x, 'y': y}
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
            ],
            [
                {'x': x, 'y': y}
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
            ]
        ]
        self.set_tolerance(10.0)

    def test_zero(self):
        """Make sure that all-zero input produces an identically zero PSF."""

        for grid, test_points in zip(self._grids, self._test_points):
            shape = (grid['y'].size, grid['x'].size)
            psf = PiecewiseBicubicPSF(
                psf_parameters={'values': numpy.zeros(shape),
                                'd_dx': numpy.zeros(shape),
                                'd_dy': numpy.zeros(shape),
                                'd2_dxdy': numpy.zeros(shape)},
                boundaries=grid
            )

            for point in test_points:
                self.assertEqual(psf(**point), 0.0)

            for point1 in test_points:
                for point2 in test_points:
                    self.assertEqual(
                        psf.integrate(left=point1['x'],
                                      bottom=point1['y'],
                                      width=point2['x'] - point1['x'],
                                      height=point2['y'] - point1['y']),
                        0.0
                    )

    def test_one(self):
        """Make sure that a PSF = 1 everywhere within the grid works."""

        for grid, test_points in zip(self._grids, self._test_points):
            shape = (grid['y'].size, grid['x'].size)
            psf = PiecewiseBicubicPSF(
                psf_parameters={'values': numpy.ones(shape),
                                'd_dx': numpy.zeros(shape),
                                'd_dy': numpy.zeros(shape),
                                'd2_dxdy': numpy.zeros(shape)},
                boundaries=grid
            )

            for point in test_points:
                self.assertEqual(psf(**point),
                                 1.0 if point_in_grid(point, grid) else 0.0,
                                 f'1.0({point["x"]:f}, {point["y"]:f})')

            for point1 in test_points:
                for point2 in test_points:
                    width = (
                        max(grid['x'][0],
                            min(grid['x'][-1], point2['x']))
                        -
                        max(grid['x'][0],
                            min(grid['x'][-1], point1['x']))
                    )
                    height = (
                        max(grid['y'][0],
                            min(grid['y'][-1], point2['y']))
                        -
                        max(grid['y'][0],
                            min(grid['y'][-1], point1['y']))
                    )
                    self.assertEqual(
                        psf.integrate(left=point1['x'],
                                      bottom=point1['y'],
                                      width=point2['x'] - point1['x'],
                                      height=point2['y'] - point1['y']),
                        width * height,
                        f"Int(1, {point1['x']:f} < x < {point2['x']:f}, "
                        f"{point1['y']:f} < y < {point2['y']:f})"
                    )


    def test_pi(self):
        """Make sure that a PSF = 1 everywhere within the grid works."""

        for grid, test_points in zip(self._grids, self._test_points):
            shape = (grid['y'].size, grid['x'].size)
            psf = PiecewiseBicubicPSF(
                psf_parameters={'values': numpy.full(shape, numpy.pi),
                                'd_dx': numpy.zeros(shape),
                                'd_dy': numpy.zeros(shape),
                                'd2_dxdy': numpy.zeros(shape)},
                boundaries=grid
            )

            for point in test_points:
                self.assertEqual(
                    psf(**point),
                    numpy.pi if point_in_grid(point, grid) else 0.0,
                    f"pi({point['x']:f}, {point['y']:f})"
                )

            for point1 in test_points:
                for point2 in test_points:
                    width = (
                        max(grid['x'][0],
                            min(grid['x'][-1], point2['x']))
                        -
                        max(grid['x'][0],
                            min(grid['x'][-1], point1['x']))
                    )
                    height = (
                        max(grid['y'][0],
                            min(grid['y'][-1], point2['y']))
                        -
                        max(grid['y'][0],
                            min(grid['y'][-1], point1['y']))
                    )
                    self.assertApprox(
                        psf.integrate(left=point1['x'],
                                      bottom=point1['y'],
                                      width=point2['x'] - point1['x'],
                                      height=point2['y'] - point1['y']),
                        width * height * numpy.pi,
                        f"Int(pi, {point1['x']:f} < x < {point2['x']:f}, "
                        f"{point1['y']:f} < y < {point2['y']:f}"
                    )

    def test_linear(self):
        """Test PSFs that is are linear functions of x and/or y work."""

        def test_evaluation(psf, slope, grid, test_points):
            """
            Test PSF evalutaion.

            Args:
                psf:    The PSF to test.

                slope:    Dictionary giving the expected coefficient of `x` or
                    `y` of the PSF (keys - `x` and `y` respectively).

                grid:    The grid on which the PSF was constructed (see
                    boundaries argument of PicewiseBicubicPSF.__init__

                test_points:    A collection of points at which to try to
                    evaluate the PSF. It is perfectly valid for the points to be
                    outside the range of the PSF.

            Returns:
                None
            """

            for point in test_points:
                if(
                        point['x'] < grid['x'][0]
                        or
                        point['y'] < grid['y'][0]
                        or
                        point['x'] > grid['x'][-1]
                        or
                        point['y'] > grid['y'][-1]
                ):
                    answer = 0.0
                else:
                    answer = slope['x'] * point['x'] + slope['y'] * point['y']
                self.assertEqual(psf(**point),
                                 answer,
                                 f"{slope['x']:f} * {point['x']:f} + "
                                 f"{slope['y']:f} * {point['y']:f}")

        def test_integration(psf, slope, grid, test_points):
            """
            Test that integrating the PSF produces the expected results.

            Args:
                See test_evaluation.

            Returns:
                None
            """

            for point1 in test_points:
                for point2 in test_points:
                    x_1 = max(
                        grid['x'][0],
                        min(grid['x'][-1], point1['x'])
                    )
                    x_2 = max(
                        grid['x'][0],
                        min(grid['x'][-1], point2['x'])
                    )
                    y_1 = max(
                        grid['y'][0],
                        min(grid['y'][-1], point1['y'])
                    )
                    y_2 = max(
                        grid['y'][0],
                        min(grid['y'][-1], point2['y'])
                    )
                    width = point2['x'] - point1['x']
                    height = point2['y'] - point1['y']
                    answer = (
                        (x_2**2 - x_1**2) * slope['x'] * (y_2 - y_1)
                        +
                        (y_2**2 - y_1**2) * slope['y'] * (x_2 - x_1)
                    ) / 2.0
                    self.assertApprox(
                        psf.integrate(left=point1['x'],
                                      bottom=point1['y'],
                                      width=width,
                                      height=height),
                        answer,
                        f"int(psf = {slope['x']:f} * x + {slope['y']:f} * y, "
                        f"{point1['x']:f} < x < {point2['x']:f}, "
                        f"{point1['y']:f} < y < {point2['y']:f}"
                    )

        slopes = [0.0, 1.0, numpy.pi]
        for grid, test_points in zip(self._grids, self._test_points):
            shape = (grid['y'].size,
                     grid['x'].size)
            for x_slope in slopes:
                for y_slope in slopes:
                    grid_x, grid_y = numpy.meshgrid(grid['x'],
                                                    grid['y'])
                    psf = PiecewiseBicubicPSF(
                        psf_parameters={
                            'values': grid_x * x_slope + grid_y * y_slope,
                            'd_dx': numpy.full(shape, x_slope),
                            'd_dy': numpy.full(shape, y_slope),
                            'd2_dxdy': numpy.zeros(shape)
                        },
                        boundaries=grid
                    )
                    psf_slope = {'x': x_slope, 'y': y_slope}

                    test_evaluation(psf, psf_slope, grid, test_points)
                    test_integration(psf, psf_slope, grid, test_points)



    #Broke this up into as many pieces as was reasonable.
    #pylint: disable=too-many-locals
    def test_random_single_patch(self):
        """Test PSFs equal to random bi-cubic polynomials."""

        def test_evaluation(psf, coef, grid, test_points):
            """
            Test PSF evalutaion.

            Args:
                psf:    The PSF to test.

                coef:    The coefficients of the PSF polynomial.

                grid:    The grid on which the PSF was constructed (see
                    boundaries argument of PicewiseBicubicPSF.__init__

                test_points:    A collection of points at which to try to
                    evaluate the PSF. It is perfectly valid for the points to be
                    outside the range of the PSF.

            Returns:
                None
            """

            for point in test_points:
                expected = 0.0
                if point_in_grid(point, grid):
                    y_term = 1.0
                    for y_pow in range(4):
                        x_term = 1.0
                        for x_pow in range(4):
                            expected += coef[y_pow, x_pow] * x_term * y_term
                            x_term *= point['x']
                        y_term *= point['y']
                got = psf(**point)
                self.assertApprox(
                    got,
                    expected,
                    f'PSF(random coef)({point["x"]:f}, {point["y"]:f})'
                )

        def test_integration(psf, coef, grid, test_points):
            """
            Test the PSF integrals over rectangular areas.

            Args:
                See test_evaluation.

            Returns:
                None
            """

            for point1 in test_points:
                for point2 in test_points:
                    expected = 0.0

                    x_1 = max(
                        grid['x'][0],
                        min(grid['x'][-1], point1['x'])
                    )
                    x_2 = max(
                        grid['x'][0],
                        min(grid['x'][-1], point2['x'])
                    )
                    y_1 = max(
                        grid['y'][0],
                        min(grid['y'][-1], point1['y'])
                    )
                    y_2 = max(
                        grid['y'][0],
                        min(grid['y'][-1], point2['y'])
                    )

                    y1_term = y_1
                    y2_term = y_2
                    for y_pow in range(1, 5):
                        x1_term = x_1
                        x2_term = x_2
                        for x_pow in range(1, 5):
                            expected += (coef[y_pow - 1, x_pow - 1]
                                         *
                                         (x2_term - x1_term) / x_pow
                                         *
                                         (y2_term - y1_term) / y_pow)
                            x1_term *= x_1
                            x2_term *= x_2
                        y1_term *= y_1
                        y2_term *= y_2

                    got = psf.integrate(left=point1['x'],
                                        bottom=point1['y'],
                                        width=point2['x'] - point1['x'],
                                        height=point2['y'] - point1['y'])

                    self.assertApprox(
                        got,
                        expected,
                        "int(PSF(random coef), "
                        f"{point1['x']:f} < x < {point2['x']:f}, "
                        f"{point1['y']:f} < y < {point2['y']:f})"
                    )


        coef = numpy.random.rand(4, 4) * numpy.pi
        for grid, test_points in zip(self._grids, self._test_points):
            print('Testing grid: ' + repr(grid))
            shape = (grid['y'].size,
                     grid['x'].size)
            grid_x, grid_y = numpy.meshgrid(grid['x'], grid['y'])
            psf_parameters = {'values': numpy.zeros(shape),
                              'd_dx': numpy.zeros(shape),
                              'd_dy': numpy.zeros(shape),
                              'd2_dxdy': numpy.zeros(shape)}
            y_term = numpy.ones(shape)
            dy_term = numpy.ones(shape)
            for y_pow in range(4):
                x_term = numpy.ones(shape)
                dx_term = numpy.ones(shape)
                for x_pow in range(4):
                    psf_parameters['values'] += (
                        coef[y_pow, x_pow] * x_term * y_term
                    )
                    psf_parameters['d_dx'] += (
                        x_pow * coef[y_pow, x_pow] * dx_term * y_term
                    )
                    psf_parameters['d_dy'] += (
                        y_pow * coef[y_pow, x_pow] * x_term * dy_term
                    )
                    psf_parameters['d2_dxdy'] += (
                        x_pow * y_pow * coef[y_pow, x_pow] * dx_term * dy_term
                    )
                    x_term *= grid_x
                    if x_pow != 0:
                        dx_term *= grid_x
                y_term *= grid_y
                if y_pow != 0:
                    dy_term *= grid_y

            psf = PiecewiseBicubicPSF(
                psf_parameters=psf_parameters,
                boundaries=grid
            )

            test_evaluation(psf, coef, grid, test_points)
            test_integration(psf, coef, grid, test_points)
    #pylint: enable=too-many-locals

if __name__ == '__main__':
    unittest.main(failfast=True)
