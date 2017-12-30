#!/usr/bin/env python3

import sys
sys.path.append('..')

from image import Image
from PiecewiseBicubicPSF import PiecewiseBicubicPSF
from Utilities import FloatTestCase

from math import floor, ceil
import unittest
import numpy

class TestImage(FloatTestCase) :
    """Make sure the Image class functions as expected."""

    def assertExpectedImage(self, got, expected, message) :
        """Check if an image is the same as expected."""

        self.assertTrue(
            numpy.vectorize(self.approx_equal, otypes = [numpy.float])(got, expected).all(),
            message
            +
            '\nGot:\n' + repr(got)
            +
            '\nExpected:\n' + repr(expected)
            +
            '\nDifference:\n' + repr(got - expected)
        )

    def setUp(self) :
        """Set tolerance to 100 x epsilon."""

        self.tolerance = 100.0

    def test_no_source(self) :
        """Test behavior of images with no sources."""

        image = Image(10, 4)
        self.assertExpectedImage(image,
                                 numpy.zeros((4, 10)),
                                 'All zero image')

        image = Image(4, 10, numpy.pi)
        self.assertExpectedImage(image,
                                 numpy.full((10, 4), numpy.pi),
                                 'All pi image.')
    
    def test_pix_aligned_source(self) :
        """Test behavior with a single source covering exact pixels."""

        images = [Image(10, 8),
                  Image(6, 8),
                  Image(6, 4),
                  Image(6, 2),
                  Image(3, 8),
                  Image(3, 4),
                  Image(3, 2)]
        shape = (2, 2)
        x0, y0 = 5.0, 3.0
        left, right = 2.0, 3.0
        down, up = 1.0, 2.0

        expected_off = 50
        expected = numpy.zeros((2 * expected_off, 2 * expected_off))
        expected[
            int(y0 - down + expected_off) : int(y0 + up + expected_off),
            int(x0 - left + expected_off) : int(x0 + right + expected_off)
        ] = 1.0

        psf = PiecewiseBicubicPSF(values = numpy.ones(shape),
                                  d_dx = numpy.zeros(shape),
                                  d_dy = numpy.zeros(shape),
                                  d2_dxdy = numpy.zeros(shape),
                                  x_boundaries = [-left, right],
                                  y_boundaries = [-down, up])

        for x_off, y_off in [(0, 0),
                             (0, -2),
                             (0, -4),
                             (0, -10),
                             (-4, 0),
                             (-4, -2),
                             (-4, -4),
                             (-4, -10),
                             (-6, 0),
                             (-6, -2),
                             (-6, -4),
                             (-6, -10),
                             (-10, 0),
                             (-10, -2),
                             (-10, -4),
                             (-10, -10),
                             (0, 2),
                             (0, 4),
                             (0, 10),
                             (4, 0),
                             (4, 2),
                             (4, 4),
                             (4, 10),
                             (6, 0),
                             (6, 2),
                             (6, 4),
                             (6, 10),
                             (10, 0),
                             (10, 2),
                             (10, 4),
                             (10, 10)] :
            for img in images :
                for subpix_map in [numpy.ones((1,1)), numpy.ones((2,3))] :
                    img.fill(0.0)
                    img.add_source(x0 + x_off,
                                   y0 + y_off,
                                   1.0,
                                   psf = psf,
                                   subpix_map = subpix_map)
                    expected_piece = expected[
                        expected_off - y_off : expected_off - y_off + img.shape[0],
                        expected_off - x_off : expected_off - x_off + img.shape[1]
                    ]
                    self.assertExpectedImage(
                        img,
                        expected_piece,
                        'Source(psf = 1, %f < x < %f, %f < y < %f)'
                        %
                        (x0 + x_off - left,
                         x0 + x_off + right,
                         y0 + y_off - down,
                         y0 + y_off + up)
                    )

    def test_misaligned_source(self) :
        """Test behavior with a single source not alinged with pixels."""

        images = [Image(10, 8),
                  Image(6, 8),
                  Image(6, 4),
                  Image(6, 2),
                  Image(3, 8),
                  Image(3, 4),
                  Image(3, 2)]

        shape = (2, 2)
        x0, y0 = 4.3, numpy.pi
        left, right = numpy.pi, numpy.e
        down, up = numpy.e / 2.0, numpy.e
        amplitude = numpy.log(10.0)
        psf = PiecewiseBicubicPSF(values = numpy.ones(shape),
                                  d_dx = numpy.zeros(shape),
                                  d_dy = numpy.zeros(shape),
                                  d2_dxdy = numpy.zeros(shape),
                                  x_boundaries = [-left, right],
                                  y_boundaries = [-down, up])
        exp_off = 50

        expected = numpy.zeros((2 * exp_off, 2 * exp_off))

        expected[
            ceil(y0 - down) + exp_off : floor(y0 + up) + exp_off,
            ceil(x0 - left) + exp_off : floor(x0 + right) + exp_off
        ] = amplitude


        expected[
            floor(y0 - down) + exp_off,
            ceil(x0 - left) + exp_off : floor(x0 + right) + exp_off
        ] = amplitude * (1.0 - (y0 - down) % 1)
        expected[
            floor(y0 + up) + exp_off,
            ceil(x0 - left) + exp_off : floor(x0 + right) + exp_off
        ] = amplitude * ((y0 + up) % 1)
        expected[
            ceil(y0 - down) + exp_off : floor(y0 + up) + exp_off,
            floor(x0 - left) + exp_off
        ] = amplitude * (1.0 - (x0 - left) % 1)
        expected[
            ceil(y0 - down) + exp_off : floor(y0 + up) + exp_off,
            floor(x0 + right) + exp_off
        ] = amplitude * ((x0 + right) % 1)


        expected[
            floor(y0 - down) + exp_off,
            floor(x0 - left) + exp_off
        ] = amplitude * (
            (1.0 - (x0 - left) % 1) * (1.0 - (y0 - down) % 1)
        )
        expected[
            floor(y0 - down) + exp_off,
            floor(x0 + right) + exp_off
        ] = amplitude * (
            ((x0 + right) % 1) * (1.0 - (y0 - down) % 1)
        )
        expected[
            floor(y0 + up) + exp_off, 
            floor(x0 - left) + exp_off
        ] = amplitude * (
            (1.0 - (x0 - left) % 1) * ((y0 + up) % 1)
        )
        expected[
            floor(y0 + up) + exp_off,
            floor(x0 + right) + exp_off
        ] = amplitude * (
            ((x0 + right) % 1) * ((y0 + up) % 1)
        )

        for x_off, y_off in [(0, 0),
                             (0, -2),
                             (0, -4),
                             (0, -10),
                             (-4, 0),
                             (-4, -2),
                             (-4, -4),
                             (-4, -10),
                             (-6, 0),
                             (-6, -2),
                             (-6, -4),
                             (-6, -10),
                             (-10, 0),
                             (-10, -2),
                             (-10, -4),
                             (-10, -10),
                             (0, 0),
                             (0, 2),
                             (0, 4),
                             (0, 10),
                             (4, 0),
                             (4, 2),
                             (4, 4),
                             (4, 10),
                             (6, 0),
                             (6, 2),
                             (6, 4),
                             (6, 10),
                             (10, 0),
                             (10, 2),
                             (10, 4),
                             (10, 10)] :
            for img in images :
                img.fill(0.0)
                img.add_source(x0 + x_off, y0 + y_off, amplitude, psf = psf)

                self.assertExpectedImage(
                    img,
                    expected[
                        exp_off - y_off : exp_off + img.shape[0] - y_off,
                        exp_off - x_off : exp_off + img.shape[1] - x_off
                    ],
                    'Source(psf = 1, %f < x < %f, %f < y %f)'
                    %
                    (x0 - left, x0 + right, y0 - down, y0 + up)
                )

if __name__ == '__main__' :
    unittest.main(failfast = True)
