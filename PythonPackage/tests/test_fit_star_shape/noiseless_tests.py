#!/usr/bin/env python3

"""Test SuperPhot's fit_star_shape module."""

import sys
import os.path
import unittest
from ctypes import c_ubyte
import numpy

_module_path = os.path.abspath(os.path.dirname(__file__))

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            _module_path,
            '..',
            '..'
        )
    )
)

#Needs to be after os.path and sys to allow adding the seach path.
#pylint: disable=wrong-import-position
from superphot import FitStarShape, BackgroundExtractor
from superphot.fake_image.piecewise_bicubic_psf import PiecewiseBicubicPSF

from tests.utilities import FloatTestCase
from tests.test_fit_star_shape.utils import make_image_and_source_list
#pylint: enable=wrong-import-position

class TestFitStarShapeNoiseless(FloatTestCase):
    """Test piecewise bicubic PSF fitting on noiseless images."""

    def check_results(self, result_tree, sources, extra_variables):
        """
        Assert that fitted PSF map and source fluxes match expectations.

        Args:
            result_tree:    The result tree containing the PSF fitting
                configuration and results.

            sources:    The sources argument used to generate the image that was
                fit. See same name argument of run_test.

            extra_variables:    A list of the names of any variables in addition
                to `x` and `y` which participate in the PSF fit.

        Returns:
            None
        """

        if 'enabled' in extra_variables:
            enabled_sources = numpy.array([src['enabled'] for src in sources],
                                          dtype=bool)
        else:
            enabled_sources = numpy.full(len(sources), True, dtype=bool)

        variables = {
            var: val[enabled_sources]
            for var, val in zip(['x', 'y'] + extra_variables,
                                result_tree.get('psffit.variables',
                                                shape=(len(sources),)))
        }

        psffit_terms = result_tree.get('psffit.terms', str)
        assert psffit_terms[0] == '{'
        assert psffit_terms[-1] == '}'

        num_sources = variables['x'].size
        num_x_boundaries = len(sources[0]['psf_args']['boundaries']['x']) - 2
        num_y_boundaries = len(sources[0]['psf_args']['boundaries']['y']) - 2

        #Using eval here is perfectly reasonable.
        #pylint: disable=eval-used
        term_list = [eval(term, variables)
                     for term in psffit_terms[1 : -1].split(',')]
        #pylint: enable=eval-used
        for term_index, term in enumerate(term_list):
            if isinstance(term, (float, int)):
                term_list[term_index] = numpy.full(num_sources, float(term))

        term_list = numpy.dstack(term_list)[0]
        coefficients = result_tree.get(
            'psffit.psfmap',
            shepe=(4,
                   len(sources[0][0]['psf_args']['boundaries']['x']) - 1,
                   len(sources[0][0]['psf_args']['boundaries']['y']) - 1,
                   len(term_list))
        )

        #Indices are: source index, variable, y boundary ind, x boundary ind
        fit_params = numpy.tensordot(term_list, coefficients, [1, 3])
        self.assertEqual(
            fit_params.shape,
            (num_sources, 4, num_x_boundaries, num_y_boundaries)
        )
        fluxes = result_tree.get('psffit.flux',
                                 shape=(len(variables['x']),))

        assert len(sources) == len(fluxes)
        for src_ind, src in enumerate(sources):
            if 'enabled' in extra_variables and not src['enabled']:
                continue

            if 'flux_backup' in src and src['flux_backup'] is not None:
                self.assertEqual(fluxes[src_ind], 0)
                fluxes[src_ind] = src['flux_backup']
            else:
                self.assertNotEqual(fluxes[src_ind], 0)

        print('fluxes before = ' + repr(fluxes))
        fluxes = fluxes[enabled_sources]
        print('fluxes after = ' + repr(fluxes))

        fit_params *= fluxes[:, numpy.newaxis, numpy.newaxis, numpy.newaxis]

        expected_params = numpy.empty(fit_params.shape)
        for src_ind, src in enumerate(sources):
            if 'enabled' in extra_variables and not src['enabled']:
                continue
            for var_ind, var_name in enumerate(['values',
                                                'd_dx',
                                                'd_dy',
                                                'd2_dxdy']):
                expected_params[src_ind, var_ind, :, :] = (
                    src['psf_args']['psf_parameters'][var_name][1:-1, 1:-1]
                )

        plus = (expected_params + fit_params)
        minus = (expected_params - fit_params)
        self.assertLess((minus * minus).sum() / (plus * plus).sum(),
                        1e-13 * minus.size,
                        msg=('Expected: ' + repr(expected_params)
                             +
                             '\n'
                             +
                             'Got: ' + repr(fit_params)))

    def run_test(self,
                 sources,
                 psffit_terms,
                 extra_variables=None):
        """
        Assert that a fit of a series of images works exactly.

        Args:
            sources:    A list of lists of dictionaries specifying the list of
                sources to fit. Each list of dictionaries specifies the sources
                to drop on a single image. Each source must contain the
                following:

                    * x:    The x coordinate of the source center.

                    * y:    The y coordinate of the source center.

                    * psf_args:    The arguments with which to create the
                      PiecewiseBicubicPSF for the source. See
                      PiecewiseBicubicPSF.__init__ for details.

            psffit_terms:    The terms on which PSF parameters depend on. See
                --psf.terms argument of the fitpsf command.

            extra_variables:    A list of the variables in addition to x and y
                that participate in the fitting terms.

        Returns:
            None
        """

        if extra_variables is None:
            extra_variables = []
        for subpixmap in [numpy.ones((1, 1))]:#,
#                           numpy.ones((2, 2)),
#                           numpy.array([[1.99, 0.01], [0.01, 1.99]]),
#                           numpy.array([[0.5, 0.5], [0.5, 2.5]]),
#                           numpy.array([[1.9], [0.1]]),
#                           numpy.array([[2.0, 0.0], [0.0, 2.0]]),
#                           numpy.array([[0.0, 0.0], [0.0, 4.0]])]:
            fit_images_and_sources = []


            measure_backgrounds = []
            for image_sources in sources:
                print('Image sources:\n' + repr(image_sources))
                image, source_list = make_image_and_source_list(
                    sources=[dict(x=src['x'],
                                  y=src['y'],
                                  psf=PiecewiseBicubicPSF(**src['psf_args']),
                                  **{var: src[var] for var in extra_variables})
                             for src in image_sources],
                    extra_variables=extra_variables,
                    subpix_map=subpixmap,
                )
                fit_images_and_sources.append(
                    (
                        image,
                        image**0.5,
                        numpy.zeros(image.shape, dtype=c_ubyte),
                        source_list
                    )
                )
                measure_backgrounds.append(
                    BackgroundExtractor(
                        image,
                        6.0,
                        7.0
                    )
                )
                measure_backgrounds[-1](
                    numpy.array([src['x'] for src in image_sources]),
                    numpy.array([src['y'] for src in image_sources])
                )

            print('Fitting for the PSF.')
            fit_star_shape = FitStarShape(
                mode='PSF',
                shape_terms=psffit_terms,
                grid=[sources[0][0]['psf_args']['boundaries']['x'],
                      sources[0][0]['psf_args']['boundaries']['y']],
                initial_aperture=5.0,
                subpixmap=subpixmap,
                smoothing=-100.0,
                max_chi2=100.0,
                pixel_rejection_threshold=100.0,
                max_abs_amplitude_change=0.0,
                max_rel_amplitude_change=1e-8,
                min_convergence_rate=-10.0,
                max_iterations=10000
            )
            result_tree = fit_star_shape.fit(
                fit_images_and_sources,
                measure_backgrounds
            )

            self.check_results(
                result_tree,
                sources,
                extra_variables
            )

    def test_single_source(self):
        """Test fitting a single source in the center of the image."""

        values = numpy.zeros((3, 3))
        d_dx = numpy.zeros((3, 3))
        d_dy = numpy.zeros((3, 3))
        d2_dxdy = numpy.zeros((3, 3))
        values[1, 1] = 1.0

        self.run_test(
            sources=[[
                dict(
                    x=15.0,
                    y=15.0,
                    psf_args=dict(
                        psf_parameters=dict(
                            values=values,
                            d_dx=d_dx,
                            d_dy=d_dy,
                            d2_dxdy=d2_dxdy
                        ),
                        boundaries=dict(x=numpy.array([-2.0, 0.0, 2.0]),
                                        y=numpy.array([-1.0, 0.0, 1.0]))
                    )
                )
            ]],
            psffit_terms='{1}'
        )

if __name__ == '__main__':
    unittest.main(failfast=True)
