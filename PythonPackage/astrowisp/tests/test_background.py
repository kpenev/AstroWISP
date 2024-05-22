#!/usr/bin/env python3

"""Test AstroWISP's background extraction."""

import unittest
import numpy

from astrowisp.tests.utilities import FloatTestCase
from astrowisp import BackgroundExtractor

class TestAnnulusBackground(FloatTestCase):
    """Test background extraction based on annuli around sources."""

    @staticmethod
    def add_flux_around_sources(source_x, source_y, radius, extra_flux, image):
        """
        Add extra flux to image pixels within a radius around sources.

        Args:
            source_x:    The x coordinates of the sources to add flux around.

            source_x:    The y coordinates of the sources to add flux around.

            radius:    Pixes with centers within this radius get extra flux.

            extra_flux:    The amount of extra flux per pixel to add.

            image:    The image to add flux to.

        Returns:
            None
        """

        for this_x, this_y in zip(source_x, source_y):
            min_x = max(0, int(numpy.ceil(this_x - radius - 0.5)))
            max_x = min(image.shape[1], int(numpy.floor(this_x + radius + 0.5)))

            min_y = max(0, int(numpy.ceil(this_y - radius - 0.5)))
            max_y = min(image.shape[0], int(numpy.floor(this_y + radius + 0.5)))

            for pixel_x in range(min_x, max_x):
                for pixel_y in range(min_y, max_y):
                    if(
                            (pixel_x + 0.5 - this_x)**2
                            +
                            (pixel_y + 0.5 - this_y)**2
                            <
                            radius**2
                    ):
                        image[pixel_y, pixel_x] += extra_flux


    def test_constant_image(self):
        """Test background extraction on an image with all pixels the same."""

        inner_radius, outer_radius = 1.99, 4.01
        for expected_value in [1.0, 10.0, 0.0]:
            for src_x, src_y in [
                    (numpy.array([5.0]), numpy.array([5.0])),
                    numpy.dstack(
                        numpy.meshgrid(
                            [2.5, 5.0, 7.5],
                            [2.5, 5.0, 7.5]
                        )
                    ).reshape(9, 2).transpose()
            ]:
                image = numpy.full(shape=(10, 10), fill_value=expected_value)
                for extra_flux in [0.0, 10.0, numpy.nan]:
                    self.add_flux_around_sources(
                        src_x,
                        src_y,
                        inner_radius,
                        extra_flux,
                        image
                    )
                    message = (
                        f'Sources at: x = {src_x!r}, y = {src_y!r}, '
                        f'extra flux = {extra_flux!r}'
                    )

                    measure_background = BackgroundExtractor(
                        image,
                        inner_radius=inner_radius,
                        outer_radius=outer_radius
                    )

                    for (
                            source_index,
                            (extracted_value, extracted_error)
                    ) in enumerate(
                        zip(
                            *measure_background(numpy.copy(src_x),
                                                numpy.copy(src_y))[:2]
                        )
                    ):
                        if source_index == 4:
                            self.assertTrue(numpy.isnan(extracted_value),
                                            message)
                            self.assertTrue(numpy.isnan(extracted_error),
                                            message)
                        else:
                            self.assertApprox(extracted_value,
                                              expected_value,
                                              message)
                            self.assertApprox(extracted_error, 0.0, message)

    def test_crowded_image(self):
        """Tests involving sources for which no BG can be determined."""

        image = numpy.ones(shape=(10, 10))

        for src_x, src_y in [
                (numpy.array([5.0]), numpy.array([5.0])),
                numpy.dstack(
                    numpy.meshgrid([2.5, 5.0, 7.5], [2.5, 5.0, 7.5])
                ).reshape(9, 2).transpose()
        ]:
            measure_background = BackgroundExtractor(image,
                                                     inner_radius=10.0,
                                                     outer_radius=15.0)

            message = f'Sources at: x = {src_x!r}, y = {src_y!r}'

            for extracted_value, extracted_error in zip(
                    *measure_background(numpy.copy(src_x),
                                        numpy.copy(src_y))[:2]
            ):
                self.assertTrue(numpy.isnan(extracted_value), message)
                self.assertTrue(numpy.isnan(extracted_error), message)

    def test_partially_crowded_image(self):
        """A test where one source is crowded and 4 are not."""

        image = numpy.ones(shape=(10, 10))

        inner_radius, outer_radius = 2.5, (3.0**0.5 + 1.0) * 1.25

        src_x = numpy.array([5.0 - 1.25,
                             5.0 + 1.25,
                             5.0 - 1.25,
                             5.0 + 1.25,
                             5.0])
        src_y = numpy.array([5.0 - 1.25,
                             5.0 - 1.25,
                             5.0 + 1.25,
                             5.0 + 1.25,
                             5.0])
        message = 'x = %f, y = %f, extra flux = %f'
        for extra_flux in [0.0, 10.0, numpy.nan]:
            self.add_flux_around_sources(
                src_x,
                src_y,
                inner_radius,
                extra_flux,
                image
            )
            measure_background = BackgroundExtractor(
                image,
                inner_radius=inner_radius,
                outer_radius=outer_radius
            )

            extracted_values, extracted_errors = measure_background(src_x,
                                                                    src_y)[:2]

            for value, source_center in zip(extracted_values[:4],
                                            zip(src_x[:4], src_y[:4])):
                self.assertApprox(value,
                                  1.0,
                                  message % (source_center + (extra_flux,)))

            for error, source_center in zip(extracted_errors[:4],
                                            zip(src_x[:4], src_y[:4])):
                self.assertApprox(error,
                                  0.0,
                                  message % (source_center + (extra_flux,)))

            self.assertTrue(numpy.isnan(extracted_values[4]),
                            message % (src_x[-1], src_y[-1], extra_flux))
            self.assertTrue(numpy.isnan(extracted_errors[4]),
                            message % (src_x[-1], src_y[-1], extra_flux))

    def test_background_gradient(self):
        """Tests where the background has a constant non-zero gradient."""

        x_gradient, y_gradient = 0.5, 1.5

        for bg_radius in [(0.99, 2.01),
                          (1.5, 2.01),
                          (1.5, 2.5),
                          (numpy.pi/3, numpy.pi/2 + 0.5)]:

            for src_x, src_y in [
                    (
                        numpy.array([5.0]),
                        numpy.array([5.0])
                    ),
                    (
                        numpy.array([2.5, 7.5, 2.5, 7.5]),
                        numpy.array([2.5, 2.5, 7.5, 7.5]),
                    )
            ]:
                half_pix = numpy.linspace(0.5, 9.5, 10)
                image = (x_gradient * half_pix[None, :]
                         +
                         y_gradient * half_pix[:, None])
                for extra_flux in [0.0, 10.0, numpy.nan]:
                    self.add_flux_around_sources(
                        src_x,
                        src_y,
                        bg_radius[0],
                        extra_flux,
                        image
                    )

                    measure_background = BackgroundExtractor(image, *bg_radius)

                    extracted_values = measure_background(src_x,
                                                          src_y)[0]
                    for eval_x, eval_y, extracted in zip(src_x,
                                                         src_y,
                                                         extracted_values):
                        self.assertApprox(
                            extracted,
                            x_gradient * eval_x + y_gradient * eval_y,
                            f'BG radius = {bg_radius[0]:f}:{bg_radius[1]:f}, '
                            f'source ({eval_x:f}, {eval_y:f}), '
                            f'extra flux = {extra_flux:f}'
                        )

    def test_background_error(self):
        """Test constant gradient background but with known error."""

        half_pix = numpy.linspace(0.5, 9.5, 10)

        source_position = numpy.arange(0.0, 10.0, 0.3 * numpy.pi)

        for source_x in source_position:
            for source_y in source_position:
                image = (0.5 * half_pix[None, :]
                         +
                         numpy.pi / 2.0 * half_pix[:, None])


                for inner in [1.0, 1.5, 2.0, numpy.pi]:
                    self.add_flux_around_sources([source_x],
                                                 [source_y],
                                                 inner,
                                                 numpy.nan,
                                                 image)

                    pixels = image.flatten()
                    pixels = pixels[numpy.logical_not(numpy.isnan(pixels))]

                    for error_confidence in numpy.linspace(0.1, 0.9, 10):
                        (
                            extracted_value,
                            extracted_error,
                            extracted_npix
                        ) = BackgroundExtractor(
                            image=image,
                            inner_radius=inner,
                            outer_radius=image.shape[0] + image.shape[1],
                            error_confidence=error_confidence
                        )(
                            numpy.array([source_x]),
                            numpy.array([source_y])
                        )

                        extracted_error /= (numpy.pi
                                            /
                                            (2.0 * (extracted_npix - 1)))**0.5

                        message = (
                            f'BG radius = {inner:f}, source '
                            f'({source_x:f}, {source_y:f}), '
                            f'confidence = {error_confidence:f}'
                        )

                        self.assertEqual(extracted_npix, pixels.size, message)

                        self.assertGreaterEqual(
                            numpy.logical_and(
                                pixels >= extracted_value - extracted_error,
                                pixels <= extracted_value + extracted_error
                            ).sum(),
                            numpy.floor(error_confidence * pixels.size + 0.499),
                            message
                        )
                        self.assertLessEqual(
                            numpy.logical_and(
                                pixels > extracted_value - extracted_error,
                                pixels < extracted_value + extracted_error
                            ).sum(),
                            numpy.floor(error_confidence * pixels.size + 0.501),
                            message
                        )

if __name__ == '__main__':
    unittest.main(failfast=False)
