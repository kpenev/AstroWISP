#!/usr/bin/env python3

import sys
sys.path.append('..')

from PiecewiseBicubicPSF import PiecewiseBicubicPSF
from Utilities import FloatTestCase
import unittest
import numpy

def point_in_grid(point, grid) :
    """Return True iff the point is within grid."""

    return (
        point['x'] >= grid['x_boundaries'][0]
        and
        point['x'] <= grid['x_boundaries'][-1]
        and
        point['y'] >= grid['y_boundaries'][0]
        and
        point['y'] <= grid['y_boundaries'][-1]
    )

class TestPiecewiseBicubicPSF(FloatTestCase) :
    """Make sure the PiecewiseBicubicPSF class functions as expected."""

    def setUp(self) :
        """Define a set of picewise PSF grids to use during tests."""

        self._grids = [
            dict(x_boundaries = numpy.array([0.0, 1.0]),
                 y_boundaries = numpy.array([0.0, 1.0])),
            dict(x_boundaries = numpy.array([0.0, 0.5]),
                 y_boundaries = numpy.array([0.0, numpy.pi / 3.0])),
            dict(x_boundaries = numpy.array([0.0, 0.5, 1.0]),
                 y_boundaries = numpy.array([0.0, 1.0])),
        ]
        self._test_points = [
            [
                dict(x = x, y = y)
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2,
                                        numpy.pi, 10)
            ],
            [
                dict(x = x, y = y)
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
            ],
            [
                dict(x = x, y = y)
                for x in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
                for y in numpy.linspace(-numpy.pi / 2, numpy.pi, 10)
            ]
        ]
        self.set_tolerance(10.0)

    def test_zero(self) :
        """Make sure that all-zero input produces an identically zero PSF."""

        for grid, test_points in zip(self._grids, self._test_points) :
            shape = (grid['y_boundaries'].size, grid['x_boundaries'].size)
            psf = PiecewiseBicubicPSF(
                values = numpy.zeros(shape),
                d_dx = numpy.zeros(shape),
                d_dy = numpy.zeros(shape),
                d2_dxdy = numpy.zeros(shape),
                **grid
            )

            for point in test_points :
                self.assertEqual(psf(**point), 0.0)

            for p1 in test_points :
                for p2 in test_points :
                    self.assertEqual(
                        psf.integrate(left = p1['x'],
                                      bottom = p1['y'],
                                      width = p2['x'] - p1['x'],
                                      height = p2['y'] - p1['y']),
                        0.0
                    )

    def test_one(self) :
        """Make sure that a PSF = 1 everywhere within the grid works."""

        for grid, test_points in zip(self._grids, self._test_points) :
            shape = (grid['y_boundaries'].size, grid['x_boundaries'].size)
            psf = PiecewiseBicubicPSF(
                values = numpy.ones(shape),
                d_dx = numpy.zeros(shape),
                d_dy = numpy.zeros(shape),
                d2_dxdy = numpy.zeros(shape),
                **grid
            )

            for point in test_points :
                self.assertEqual(psf(**point), 
                                 1.0 if point_in_grid(point, grid) else 0.0,
                                 '1.0(%f, %f)' % (point['x'], point['y']))

            for p1 in test_points :
                for p2 in test_points :
                    width = (
                        max(grid['x_boundaries'][0],
                            min(grid['x_boundaries'][-1], p2['x']))
                        - 
                        max(grid['x_boundaries'][0],
                            min(grid['x_boundaries'][-1], p1['x']))
                    )
                    height = (
                        max(grid['y_boundaries'][0],
                            min(grid['y_boundaries'][-1], p2['y']))
                        -
                        max(grid['y_boundaries'][0],
                            min(grid['y_boundaries'][-1], p1['y']))
                    )
                    self.assertEqual(
                        psf.integrate(left = p1['x'],
                                      bottom = p1['y'],
                                      width = p2['x'] - p1['x'],
                                      height = p2['y'] - p1['y']),
                        width * height,
                        'Int(1, %f < x < %f, %f < y < %f)'
                        %
                        (p1['x'], p2['x'], p1['y'], p2['y'])
                    )


    def test_pi(self) :
        """Make sure that a PSF = 1 everywhere within the grid works."""

        for grid, test_points in zip(self._grids, self._test_points) :
            shape = (grid['y_boundaries'].size, grid['x_boundaries'].size)
            psf = PiecewiseBicubicPSF(
                values = numpy.full(shape, numpy.pi),
                d_dx = numpy.zeros(shape),
                d_dy = numpy.zeros(shape),
                d2_dxdy = numpy.zeros(shape),
                **grid
            )

            for point in test_points :
                self.assertEqual(
                    psf(**point),
                    numpy.pi if point_in_grid(point, grid) else 0.0,
                    'pi(%f, %f)' % (point['x'], point['y'])
                )

            for p1 in test_points :
                for p2 in test_points :
                    width = (
                        max(grid['x_boundaries'][0],
                            min(grid['x_boundaries'][-1], p2['x']))
                        - 
                        max(grid['x_boundaries'][0],
                            min(grid['x_boundaries'][-1], p1['x']))
                    )
                    height = (
                        max(grid['y_boundaries'][0],
                            min(grid['y_boundaries'][-1], p2['y']))
                        -
                        max(grid['y_boundaries'][0],
                            min(grid['y_boundaries'][-1], p1['y']))
                    )
                    self.assertApprox(
                        psf.integrate(left = p1['x'],
                                      bottom = p1['y'],
                                      width = p2['x'] - p1['x'],
                                      height = p2['y'] - p1['y']),
                        width * height * numpy.pi,
                        'Int(pi, %f < x < %f, %f < y < %f'
                        %
                        (p1['x'], p2['x'], p1['y'], p2['y'])
                    )

    def test_linear(self) :
        """Test PSFs that is are linear functions of x and/or y work."""

        slopes = [0.0, 1.0, numpy.pi]
        for grid, test_points in zip(self._grids, self._test_points) :
            shape = (grid['y_boundaries'].size,
                     grid['x_boundaries'].size)
            for x_slope in slopes :
                for y_slope in slopes :
                    grid_x, grid_y = numpy.meshgrid(grid['x_boundaries'],
                                                    grid['y_boundaries'])
                    psf = PiecewiseBicubicPSF(
                        values = grid_x * x_slope + grid_y * y_slope,
                        d_dx = numpy.full(shape, x_slope),
                        d_dy = numpy.full(shape, y_slope),
                        d2_dxdy = numpy.zeros(shape),
                        **grid
                    )

                    for point in test_points :
                        if(
                                point['x'] < grid['x_boundaries'][0]
                                or
                                point['y'] < grid['y_boundaries'][0]
                                or
                                point['x'] > grid['x_boundaries'][-1]
                                or
                                point['y'] > grid['y_boundaries'][-1]
                        ) : 
                            answer = 0.0
                        else :
                            answer = x_slope * point['x'] + y_slope * point['y']
                        self.assertEqual(
                            psf(**point),
                            answer,
                            '%f * %f + %f * %f'
                            %
                            (x_slope, point['x'], y_slope, point['y'])
                        )

                    for p1 in test_points :
                        for p2 in test_points :
                            x1 = max(
                                grid['x_boundaries'][0],
                                min(grid['x_boundaries'][-1], p1['x'])
                            )
                            x2 = max(
                                grid['x_boundaries'][0],
                                min(grid['x_boundaries'][-1], p2['x'])
                            )
                            y1 = max(
                                grid['y_boundaries'][0],
                                min(grid['y_boundaries'][-1], p1['y'])
                            )
                            y2 = max(
                                grid['y_boundaries'][0],
                                min(grid['y_boundaries'][-1], p2['y'])
                            )
                            width = p2['x'] - p1['x']
                            height = p2['y'] - p1['y']
                            answer = (
                                (x2**2 - x1**2) * x_slope * (y2 - y1)
                                +
                                (y2**2 - y1**2) * y_slope * (x2 - x1)
                            ) / 2.0
                            self.assertApprox(
                                psf.integrate(left = p1['x'],
                                              bottom = p1['y'],
                                              width = width,
                                              height = height),
                                answer,
                                (
                                    'int(psf = %f * x + %f * y, '
                                    '%f < x < %f, '
                                    '%f < y < %f'
                                )
                                %
                                (x_slope,
                                 y_slope,
                                 p1['x'],
                                 p2['x'],
                                 p1['y'],
                                 p2['y'])
                            )

    def test_random_single_patch(self) :
        """Test PSFs equal to random bi-cubic polynomials."""

        coef = numpy.random.rand(4, 4) * numpy.pi
        for grid, test_points in zip(self._grids, self._test_points) :
            print('Testing grid: ' + repr(grid))
            shape = (grid['y_boundaries'].size,
                     grid['x_boundaries'].size)
            grid_x, grid_y = numpy.meshgrid(grid['x_boundaries'],
                                            grid['y_boundaries'])
            values = numpy.zeros(shape)
            d_dx = numpy.zeros(shape)
            d_dy = numpy.zeros(shape)
            d2_dxdy = numpy.zeros(shape)
            y_term = numpy.ones(shape)
            dy_term = numpy.ones(shape)
            for y_pow in range(4) :
                x_term = numpy.ones(shape)
                dx_term = numpy.ones(shape)
                for x_pow in range(4) :
                    values += coef[y_pow, x_pow] * x_term * y_term
                    d_dx += (x_pow * coef[y_pow, x_pow]
                             *
                             dx_term
                             *
                             y_term)
                    d_dy += (y_pow * coef[y_pow, x_pow]
                             *
                             x_term
                             *
                             dy_term)
                    d2_dxdy += (x_pow * y_pow * coef[y_pow, x_pow]
                                *
                                dx_term
                                *
                                dy_term)
                    x_term *= grid_x
                    if x_pow != 0 : dx_term *= grid_x
                y_term *= grid_y
                if y_pow != 0 : dy_term *= grid_y

            psf = PiecewiseBicubicPSF(
                values = values,
                d_dx = d_dx,
                d_dy = d_dy,
                d2_dxdy = d2_dxdy,
                **grid
            )

            for point in test_points :
                expected = 0.0
                if point_in_grid(point, grid) :
                    y_term = 1.0
                    for y_pow in range(4) :
                        x_term = 1.0
                        for x_pow in range(4) :
                            expected += coef[y_pow, x_pow] * x_term * y_term
                            x_term *= point['x']
                        y_term *= point['y']
                got = psf(**point)
                self.assertApprox(
                    got,
                    expected,
                    'PSF(random coef)(%(x)f, %(y)f)' % point
                )

            for p1 in test_points :
                for p2 in test_points :
                    expected = 0.0

                    x1 = max(
                        grid['x_boundaries'][0],
                        min(grid['x_boundaries'][-1], p1['x'])
                    )
                    x2 = max(
                        grid['x_boundaries'][0],
                        min(grid['x_boundaries'][-1], p2['x'])
                    )
                    y1 = max(
                        grid['y_boundaries'][0],
                        min(grid['y_boundaries'][-1], p1['y'])
                    )
                    y2 = max(
                        grid['y_boundaries'][0],
                        min(grid['y_boundaries'][-1], p2['y'])
                    )

                    y1_term = y1
                    y2_term = y2
                    for y_pow in range(1, 5) :
                        x1_term = x1
                        x2_term = x2
                        for x_pow in range(1, 5) :
                            expected += (coef[y_pow - 1 , x_pow - 1]
                                         *
                                         (x2_term - x1_term) / x_pow
                                         *
                                         (y2_term - y1_term) / y_pow)
                            x1_term *= x1
                            x2_term *= x2
                        y1_term *= y1
                        y2_term *= y2

                    got = psf.integrate(left = p1['x'],
                                        bottom = p1['y'],
                                        width = p2['x'] - p1['x'],
                                        height = p2['y'] - p1['y'])

                    self.assertApprox(
                        got,
                        expected,
                        'int(PSF(random coef), %f < x < %f, %f < y < %f)'
                        %
                        (p1['x'], p2['x'], p1['y'], p2['y'])
                    )

if __name__ == '__main__' :
    unittest.main(failfast = True)
