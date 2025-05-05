#!/usr/bin/env python3

"""Unittest test case for the astrowisp.fake_image.image.Image class."""

#Similar imports triggers pylint complaint, disabling.
from math import floor, ceil

import unittest
import numpy

from astrowisp.tests.utilities import FloatTestCase
from astrowisp.fake_image.image import Image
from astrowisp.fake_image.piecewise_bicubic_psf import PiecewiseBicubicPSF

class TestImage(FloatTestCase):
    """Make sure the Image class functions as expected."""

    #Follow the namign convention for standard unittest asserts and setUp.
    #pylint: disable=invalid-name
    def assertExpectedImage(self, got, expected, message):
        """Check if an image is the same as expected."""

        self.assertTrue(
            numpy.vectorize(
                self.approx_equal,
                otypes=[float]
            )(got, expected).all(),
            (
                message
                +
                '\nGot:\n' + repr(got)
                +
                '\nExpected:\n' + repr(expected)
                +
                '\nDifference:\n' + repr(got - expected)
            )
        )

    def setUp(self):
        """Set tolerance to 100 x epsilon."""

        self.set_tolerance(100.0)
    #pylint: enable=invalid-name

    def test_no_source(self):
        """Test behavior of images with no sources."""

        image = Image(10, 4)
        self.assertExpectedImage(image,
                                 numpy.zeros((4, 10)),
                                 'All zero image')

        image = Image(4, 10, numpy.pi)
        self.assertExpectedImage(image,
                                 numpy.full((10, 4), numpy.pi),
                                 'All pi image.')

    #The many local variables help with readability
    #pylint: disable=too-many-locals
    def test_pix_aligned_source(self):
        """Test behavior with a single source covering exact pixels."""

        images = [Image(10, 8),
                  Image(6, 8),
                  Image(6, 4),
                  Image(6, 2),
                  Image(3, 8),
                  Image(3, 4),
                  Image(3, 2)]
        shape = (2, 2)
        x_0, y_0 = 5.0, 3.0
        left, right = 2.0, 3.0
        #pylint false positive: up is a perfectly reasonable name.
        #pylint: disable=invalid-name
        down, up = 1.0, 2.0
        #pylint: enable=invalid-name

        expected_off = 50
        expected = numpy.zeros((2 * expected_off, 2 * expected_off))
        expected[
            int(y_0 - down + expected_off) : int(y_0 + up + expected_off),
            int(x_0 - left + expected_off) : int(x_0 + right + expected_off)
        ] = 1.0

        psf = PiecewiseBicubicPSF(
            psf_parameters={'values': numpy.ones(shape),
                            'd_dx': numpy.zeros(shape),
                            'd_dy': numpy.zeros(shape),
                            'd2_dxdy': numpy.zeros(shape)},
            boundaries={'x': [-left, right],
                        'y': [-down, up]}
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
                             (10, 10)]:
            for img in images:
                for subpix_map in [numpy.ones((1, 1)), numpy.ones((2, 3))]:
                    img.fill(0.0)
                    img.add_source(x_0 + x_off,
                                   y_0 + y_off,
                                   1.0,
                                   psf=psf,
                                   subpix_map=subpix_map)
                    expected_piece = expected[
                        expected_off - y_off
                        :
                        expected_off - y_off + img.shape[0]
                        ,
                        expected_off - x_off
                        :
                        expected_off - x_off + img.shape[1]
                    ]
                    self.assertExpectedImage(
                        img,
                        expected_piece,
                        f'Source(psf = 1, {x_0 + x_off - left:f}'
                        ' < x < '
                        f'{x_0 + x_off + right:f}, '
                        f'{y_0 + y_off - down:f}'
                        ' < y < '
                        f'{y_0 + y_off + up:f})'
                    )

    def test_misaligned_source(self):
        """Test behavior with a single source not alinged with pixels."""

        images = [Image(10, 8),
                  Image(6, 8),
                  Image(6, 4),
                  Image(6, 2),
                  Image(3, 8),
                  Image(3, 4),
                  Image(3, 2)]

        shape = (2, 2)
        x_0, y_0 = 4.3, numpy.pi
        left, right = numpy.pi, numpy.e
        #pylint false positive - up is a perfectly reasonable variable name
        #pylint: disable=invalid-name
        down, up = numpy.e / 2.0, numpy.e
        #pylint: enable=invalid-name
        amplitude = numpy.log(10.0)
        psf = PiecewiseBicubicPSF(
            psf_parameters={'values': numpy.ones(shape),
                            'd_dx': numpy.zeros(shape),
                            'd_dy': numpy.zeros(shape),
                            'd2_dxdy': numpy.zeros(shape)},
            boundaries={'x': [-left, right],
                        'y': [-down, up]}
        )
        exp_off = 50

        expected = numpy.zeros((2 * exp_off, 2 * exp_off))

        expected[
            ceil(y_0 - down) + exp_off : floor(y_0 + up) + exp_off,
            ceil(x_0 - left) + exp_off : floor(x_0 + right) + exp_off
        ] = amplitude


        expected[
            floor(y_0 - down) + exp_off,
            ceil(x_0 - left) + exp_off : floor(x_0 + right) + exp_off
        ] = amplitude * (1.0 - (y_0 - down) % 1)
        expected[
            floor(y_0 + up) + exp_off,
            ceil(x_0 - left) + exp_off : floor(x_0 + right) + exp_off
        ] = amplitude * ((y_0 + up) % 1)
        expected[
            ceil(y_0 - down) + exp_off : floor(y_0 + up) + exp_off,
            floor(x_0 - left) + exp_off
        ] = amplitude * (1.0 - (x_0 - left) % 1)
        expected[
            ceil(y_0 - down) + exp_off : floor(y_0 + up) + exp_off,
            floor(x_0 + right) + exp_off
        ] = amplitude * ((x_0 + right) % 1)


        expected[
            floor(y_0 - down) + exp_off,
            floor(x_0 - left) + exp_off
        ] = amplitude * (
            (1.0 - (x_0 - left) % 1) * (1.0 - (y_0 - down) % 1)
        )
        expected[
            floor(y_0 - down) + exp_off,
            floor(x_0 + right) + exp_off
        ] = amplitude * (
            ((x_0 + right) % 1) * (1.0 - (y_0 - down) % 1)
        )
        expected[
            floor(y_0 + up) + exp_off,
            floor(x_0 - left) + exp_off
        ] = amplitude * (
            (1.0 - (x_0 - left) % 1) * ((y_0 + up) % 1)
        )
        expected[
            floor(y_0 + up) + exp_off,
            floor(x_0 + right) + exp_off
        ] = amplitude * (
            ((x_0 + right) % 1) * ((y_0 + up) % 1)
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
                             (10, 10)]:
            for img in images:
                img.fill(0.0)
                img.add_source(x_0 + x_off, y_0 + y_off, amplitude, psf=psf)

                self.assertExpectedImage(
                    img,
                    expected[
                        exp_off - y_off : exp_off + img.shape[0] - y_off,
                        exp_off - x_off : exp_off + img.shape[1] - x_off
                    ],
                    f'Source(psf = 1, {x_0 - left:f} < x < {x_0 + right:f}, '
                    f'{y_0 - down:f} < y {y_0 + up:f})'
                )
    #pylint: enable=too-many-locals

if __name__ == '__main__':
    unittest.main(failfast=True)
