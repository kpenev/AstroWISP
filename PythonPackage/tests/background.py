#!/usr/bin/env python3

"""Test SuperPhot's background extraction."""

import sys
import os.path
import numpy

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            '..'
        )
    )
)

#Needs to be after os.path and sys to allow adding the seach path.
#pylint: disable=wrong-import-position
from tests.utilities import FloatTestCase
from superphot import BackgroundExtractor
#pylint: enable=wrong-import-position

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
            max_x = min(image.shape[1], int(numpy.floor(this_x + radius - 0.5)))

            min_y = max(0, int(numpy.ceil(this_y - radius - 0.5)))
            max_y = min(image.shape[0], int(numpy.floor(this_y + radius - 0.5)))

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

        measure_background = BackgroundExtractor(inner_radius=2.0,
                                                 outer_radius=4.0)
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
                        measure_background.inner_radius,
                        extra_flux,
                        image
                    )
                    message = (
                        'Sources at: x = %s, y = %s, extra flux = %s'
                        %
                        (repr(src_x), repr(src_y), extra_flux)
                    )

                    for extracted_value, extracted_error in zip(
                            *measure_background(image, src_x, src_y)
                    ):
                        self.assertApprox(extracted_value,
                                          expected_value,
                                          message)
                        self.assertApprox(extracted_error, 0.0, message)

    def test_crowded_image(self):
        """Tests involving sources for which no BG can be determined."""

        measure_background = BackgroundExtractor(inner_radius=10.0,
                                                 outer_radius=15.0)
        image = numpy.ones(shape=(10, 10))
        for src_x, src_y in [
                (numpy.array([5.0]), numpy.array([5.0])),
                numpy.dstack(
                    numpy.meshgrid([2.5, 5.0, 7.5], [2.5, 5.0, 7.5])
                ).reshape(9, 2).transpose()
        ]:
            message = ('Sources at: x = %s, y = %s'
                       %
                       (repr(src_x), repr(src_y)))

            for extracted_value, extracted_error in zip(
                    *measure_background(image, src_x, src_y)
            ):
                self.assertTrue(numpy.isnan(extracted_value), message)
                self.assertTrue(numpy.isnan(extracted_error), message)

    def test_partially_crowded_image(self):
        """A test where one source is crowded and 4 are not."""

        measure_background = BackgroundExtractor(
            inner_radius=5.0,
            outer_radius=(3.0**0.5 + 1.0) * 2.5
        )
        image = numpy.ones(shape=(10, 10))
        src_x = numpy.array([2.5, 7.5, 2.5, 7.5, 5.0])
        src_y = numpy.array([2.5, 2.5, 7.5, 7.5, 5.0])
        message = 'x = %f, y = %f, extra flux = %f'
        for extra_flux in [0.0, 10.0, numpy.nan]:
            self.add_flux_around_sources(
                src_x,
                src_y,
                measure_background.inner_radius,
                extra_flux,
                image
            )
            extracted_values, extracted_errors = measure_background(image,
                                                                    src_x,
                                                                    src_y)

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

        for bg_radius in [(1.0, 2.0),
                          (1.5, 2.0),
                          (1.5, 2.5),
                          (numpy.pi/2, 5.0 - numpy.pi/2)]:
            measure_background = BackgroundExtractor(*bg_radius)

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
                        measure_background.inner_radius,
                        extra_flux,
                        image
                    )
                    extracted_values = measure_background(image,
                                                          src_x,
                                                          src_y)[0]
                    for eval_x, eval_y, extracted in zip(src_x,
                                                         src_y,
                                                         extracted_values):
                        self.assertApprox(
                            extracted,
                            x_gradient * eval_x + y_gradient * eval_y,
                            (
                                'BG radius = %f:%f, source (%f, %f), '
                                'extra flux = %f'
                            )
                            %
                            (bg_radius + (eval_x, eval_y, extra_flux))
                        )

    def test_background_error(self):
        """Test constant gradient background but with known error."""

        half_pix = numpy.linspace(0.5, 9.5, 10)

        source_position = numpy.arange(0.0, 10.0, 0.1 * numpy.pi)

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

                    for error_confidence in numpy.linspace(0.1, 0.9, 10.0):
                        extracted_value, extracted_error = BackgroundExtractor(
                            inner_radius=inner,
                            outer_radius=image.shape[0] + image.shape[1],
                            error_confidence=error_confidence
                        )(image, numpy.array([source_x]), numpy.array(source_y))

                        min_range = extracted_value - extracted_error
                        max_range = extracted_value + extracted_error

                        num_in_range = numpy.logical_and(
                            pixels > min_range,
                            pixels < max_range
                        ).sum()

                        self.assertEqual(
                            num_in_range,
                            numpy.ceil(error_confidence * pixels.size()),
                            'BG radius = %f, source (%f, %f), confidence = %f'
                            %
                            (inner, source_x, source_y, error_confidence)
                        )
