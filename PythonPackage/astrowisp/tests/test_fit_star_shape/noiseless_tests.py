#!/usr/bin/env python3

"""Test AstroWISP's fit_star_shape module."""

import os.path
import unittest
from ctypes import c_ubyte
import numpy


from astrowisp import FitStarShape, BackgroundExtractor

from astrowisp.tests.utilities import FloatTestCase
from astrowisp.tests.test_fit_star_shape.utils import\
    make_image_and_source_list,\
    evaluate_psffit_terms

from astrowisp.fake_image.piecewise_bicubic_psf import PiecewiseBicubicPSF

class TestFitStarShapeNoiseless(FloatTestCase):
    """Test piecewise bicubic PSF fitting on noiseless images."""

    def create_debug_files(self,
                           image,
                           source_list,
                           fit_config,
                           sub_image=None):
        """
        Create the pair of files used by the C test of PSF fitting.

        Args:
            image (2D numpy array):    The image being fit.

            source_list:    The list of sources participating in the fit.

            fit_config:    The :attr:`FitStarShape.configuration` of the PSF
                fitting object used for fitting.

            sub_image:    The index of the image within the list of images being
                fit simultaneously.
        """

        fname_start = (
            os.path.expanduser('~/projects/git/AstroWISP/src/debug/')
            +
            self.id().rsplit('.', 1)[1]
            +
            '_' + str(sub_image)
        )

        with open(fname_start + '_config.txt',
                  'w',
                  encoding='utf-8') as test_config:
            for param_value in fit_config.items():
                #Needed for debugging purposes
                #pylint: disable=protected-access
                formatted_config = FitStarShape._format_config(param_value)
                #pylint: enable=protected-access
                if formatted_config:
                    test_config.write(formatted_config[0].decode()
                                      +
                                      ' = '
                                      +
                                      formatted_config[1].decode()
                                      +
                                      '\n')

        with open(fname_start + '_image.txt',
                  'w',
                  encoding='utf-8') as test_image:
            test_image.write(str(image.shape[1])
                             +
                             ' '
                             +
                             str(image.shape[0])
                             +
                             '\n')
            for value in image.flatten():
                test_image.write('\n' + repr(value))

        with open(fname_start + '_sources.txt',
                  'w',
                  encoding='utf-8') as test_sources:
            for var in source_list.dtype.names:
                test_sources.write(f'{var:25s}')
            test_sources.write('\n')
            for source in source_list:
                for var in source_list.dtype.names:
                    if var == 'ID':
                        test_sources.write(f'{source[var].decode():25s}')
                    else:
                        test_sources.write(f'{source[var]:25.16e}')
                test_sources.write('\n')

    #TODO: See if breaking up makes sense
    #pylint: disable=too-many-locals
    def check_results(self, result_tree, image_index, sources, num_terms):
        """
        Assert that fitted PSF map and source fluxes match expectations.

        Args:
            result_tree:    The result tree containing the PSF fitting
                configuration and results.

            image_index:    The index of the image for which to check results
                within the result tree (the same as the index when fitting was
                called).

            sources:    The sources argument used to generate the image that was
                fit. See same name argument of run_test.

            num_terms(int):    The number of terms the PSF map depends on.

        Returns:
            None
        """

        enabled_sources = numpy.array(
            [src.get('enabled', True) for src in sources],
            dtype=bool
        )
        num_enabled_sources = enabled_sources.sum()

        psffit_terms = result_tree.get(f'psffit.terms.{image_index:d}',
                                       shape=(len(sources), num_terms))
        psffit_terms = psffit_terms[enabled_sources]

        num_x_boundaries = len(sources[0]['psf_args']['boundaries']['x']) - 2
        num_y_boundaries = len(sources[0]['psf_args']['boundaries']['y']) - 2

        coefficients = result_tree.get(
            'psffit.psfmap',
            shape=(4,
                   len(sources[0]['psf_args']['boundaries']['x']) - 2,
                   len(sources[0]['psf_args']['boundaries']['y']) - 2,
                   num_terms)
        )

        #Indices are: source index, variable, y boundary ind, x boundary ind
        fit_params = numpy.tensordot(psffit_terms, coefficients, [1, 3])
        self.assertEqual(
            fit_params.shape,
            (num_enabled_sources, 4, num_x_boundaries, num_y_boundaries)
        )
        fluxes = result_tree.get('psffit.flux.' + str(image_index),
                                 shape=(len(sources),))

        for src_ind, src in enumerate(sources):
            if not enabled_sources[src_ind]:
                continue

            if 'flux_backup' in src and src['flux_backup'] is not None:
                self.assertEqual(fluxes[src_ind], 0)
                fluxes[src_ind] = src['flux_backup']
            else:
                self.assertNotEqual(fluxes[src_ind], 0)

        fluxes = fluxes[enabled_sources]

        fit_params *= fluxes[:, numpy.newaxis, numpy.newaxis, numpy.newaxis]

        expected_params = numpy.empty(fit_params.shape)

        enabled_ind = 0
        for src in sources:
            if not src.get('enabled', True):
                continue
            for var_ind, var_name in enumerate(['values',
                                                'd_dx',
                                                'd_dy',
                                                'd2_dxdy']):
                expected_params[enabled_ind, var_ind, :, :] = (
                    src['psf_args']['psf_parameters'][var_name][1:-1, 1:-1]
                )
            enabled_ind += 1

        plus = expected_params + fit_params
        minus = expected_params - fit_params
        self.assertLess((minus * minus).sum() / (plus * plus).sum(),
                        1e-8 * minus.size,
                        msg=('Expected: ' + repr(expected_params)
                             +
                             '\n'
                             +
                             'Got: ' + repr(fit_params)))
    #pylint: enable=too-many-locals

    def run_test(self,
                 sources,
                 psffit_terms):
        """
        Assert that a fit of a series of images works exactly.

        Args:
            sources:    A list of lists of dictionaries specifying the list of
                sources to fit. Each list of dictionaries specifies the sources
                to drop on a single image. Each source must contain at least the
                following, as well as additional variables needed to evaluate
                `psffit_terms`:

                    * x:    The x coordinate of the source center.

                    * y:    The y coordinate of the source center.

                    * psf_args:    The arguments with which to create the
                      PiecewiseBicubicPSF for the source. See
                      PiecewiseBicubicPSF.__init__ for details.

            psffit_terms[str]:    List of expressions involving entries in
                sources the PSF map will depend linearly on (e.g. `'x**2 + y'`).
                See :mod:`asteval` documentation for a list of available
                functions.

        Returns:
            None
        """

        for subpixmap in [
                numpy.ones((1, 1)),
                numpy.ones((1, 2)),
                numpy.ones((2, 1)),
                numpy.ones((2, 2)),
                numpy.array([[1.99, 0.01], [0.01, 1.99]]),
                numpy.array([[0.5, 0.5], [0.5, 2.5]]),
                numpy.array([[1.9], [0.1]]),
                numpy.array([[2.0, 0.0], [0.0, 2.0]]),
                numpy.array([[0.0, 0.0], [0.0, 4.0]])
        ]:

            #print('Fitting for the PSF.')
            fit_star_shape = FitStarShape(
                mode='PSF',
                grid=[sources[0][0]['psf_args']['boundaries']['x'],
                      sources[0][0]['psf_args']['boundaries']['y']],
                initial_aperture=5.0,
                subpixmap=subpixmap,
                smoothing=-100.0,
                max_chi2=100.0,
                pixel_rejection_threshold=100.0,
                max_abs_amplitude_change=0.0,
                max_rel_amplitude_change=1e-13,
                min_convergence_rate=-10.0,
                max_iterations=10000,
                bg_min_pix=3
            )

            fit_images_and_sources = []
            measure_backgrounds = []
            for sub_image, image_sources in enumerate(sources):
                #print(f'Sub-image #{sub_image:d} sources:\n')
                #for src in image_sources:
                #    print('\t' + repr(src) + '\n')
                psf_sources = [
                    {
                        'x': src['x'],
                        'y': src['y'],
                        'enabled': src.get('enabled', True),
                        'psf': PiecewiseBicubicPSF(**src['psf_args'])
                    }
                    for src in image_sources
                ]
                image, source_list = make_image_and_source_list(
                    sources=psf_sources,
                    subpix_map=subpixmap,
                )
                fit_images_and_sources.append(
                    (
                        image,
                        image**0.5,
                        numpy.zeros(image.shape, dtype=c_ubyte),
                        source_list,
                        evaluate_psffit_terms(image_sources, psffit_terms)
                    )
                )
                measure_backgrounds.append(
                    BackgroundExtractor(
                        fit_images_and_sources[-1][0],
                        6.0,
                        13.0
                    )
                )
                measure_backgrounds[-1](
                    numpy.array([src['x'] for src in image_sources]),
                    numpy.array([src['y'] for src in image_sources])
                )
#                self.create_debug_files(image,
#                                        source_list,
#                                        fit_star_shape.configuration,
#                                        sub_image)

            #print(80*'=')
            #print('Fitting for star shape')
            #print(80*'=')
            result_tree = fit_star_shape.fit(
                fit_images_and_sources,
                measure_backgrounds
            )

            for image_index, image_sources in enumerate(sources):
                self.check_results(
                    result_tree,
                    image_index,
                    image_sources,
                    len(psffit_terms)
                )
                #print('Finished checking results for image '
                #      +
                #      str(image_index))

    def test_single_source(self):
        """Test fitting a single source in the center of the image."""

        values = numpy.zeros((3, 3))
        d_dx = numpy.zeros((3, 3))
        d_dy = numpy.zeros((3, 3))
        d2_dxdy = numpy.zeros((3, 3))
        values[1, 1] = 1.0

        self.run_test(
            sources=[[
                {
                    'x': 15.0,
                    'y': 15.0,
                    'psf_args': {
                        'psf_parameters': {
                            'values': values,
                            'd_dx': d_dx,
                            'd_dy': d_dy,
                            'd2_dxdy': d2_dxdy
                        },
                        'boundaries': {'x': numpy.array([-2.0, 0.0, 2.0]),
                                       'y': numpy.array([-1.0, 0.0, 1.0])}
                    }
                }
            ]],
            psffit_terms=['1']
        )

    def test_isolated_sources(self):
        """Test fitting an image containing 8 well isolated sources."""

        psf_parameters = {'values': numpy.zeros((3, 3)),
                          'd_dx': numpy.zeros((3, 3)),
                          'd_dy': numpy.zeros((3, 3)),
                          'd2_dxdy': numpy.zeros((3, 3))}
        boundaries = {'x': numpy.array([-2.0, 0.0, 2.0]),
                      'y': numpy.array([-1.4, 0.0, 1.4])}

        sources = []

        psf_parameters['values'][1, 1] = 1.0
        sources.append(
            {
                'x': 15.0,
                'y': 15.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 1.0
        sources.append(
            {
                'x': 45.0,
                'y': 15.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(
            {
                'x': 15.0,
                'y': 45.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.75
        psf_parameters['d_dy'][1, 1] = 1.00
        sources.append(
            {
                'x': 37.5,
                'y': 45.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.5
        psf_parameters['d_dy'][1, 1] = 0.0
        sources.append(
            {
                'x': 30.0,
                'y': 15.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.0
        psf_parameters['d_dy'][1, 1] = 0.5
        sources.append(
            {
                'x': 15.0,
                'y': 30.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.5
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(
            {
                'x': 30.0,
                'y': 45.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 0.5
        sources.append(
            {
                'x': 45.0,
                'y': 30.0,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        self.run_test(sources=[sources],
                      psffit_terms=['1', 'x', 'y'])

    def test_two_overlapping_sources(self):
        """Test fitting an image containing 2 sources all overlapping."""

        psf_args = {
            'psf_parameters': {
                'values': numpy.zeros((3, 3)),
                'd_dx': numpy.zeros((3, 3)),
                'd_dy': numpy.zeros((3, 3)),
                'd2_dxdy': numpy.zeros((3, 3))
            },
            'boundaries': {
                'x': numpy.array([-1.0, 0.0, 1.0]),
                'y': numpy.array([-1.4, 0.0, 1.4])
            }
        }

        psf_args['psf_parameters']['values'][1, 1] = 1.0

        sources = [
            {
                'x': 14.53,
                'y': 14.02,
                'psf_args': psf_args
            },
            {
                'x': 16.51,
                'y': 14.03,
                'psf_args': psf_args
            }
        ]
        self.run_test(sources=[sources], psffit_terms=['1'])

    def test_four_overlapping_sources(self):
        """Test fitting an image containing 4 overlapping sources."""

        psf_parameters = {
            'values': numpy.zeros((3, 3)),
            'd_dx': numpy.zeros((3, 3)),
            'd_dy': numpy.zeros((3, 3)),
            'd2_dxdy': numpy.zeros((3, 3))
        }

        boundaries = {
            'x': numpy.array([-2.2, 0.0, 2.2]),
            'y': numpy.array([-1.4, 0.0, 1.4])
        }

        psf_parameters['values'][1, 1] = 1.0

        sources = [
            {
                'x': 12.5,
                'y': 13.5,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        ]


        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        sources.append(
            {
                'x': 16.5,
                'y': 13.5,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )

        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(
            {
                'x': 12.5,
                'y': 15.5,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )

        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(
            {
                'x': 16.5,
                'y': 15.5,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )

        for src in sources:
            src['enabled'] = True

        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d2_dxdy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 2.0
        psf_parameters['d2_dxdy'][1, 1] = 3.0
        sources.append(
            {
                'x': 4.0,
                'y': 4.0,
                'enabled': False,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        sources.append(
            {
                'x': 25.0,
                'y': 25.0,
                'enabled': False,
                'psf_args': {'psf_parameters': dict(psf_parameters),
                             'boundaries': boundaries}
            }
        )
        self.run_test(sources=[sources],
                      psffit_terms=['1', 'x*x', 'y*y'])

    def test_multi_image_with_extra_var(self):
        """Test fitting a series of 5 images and non-position variables."""

        boundaries = {'x': numpy.array([-3.02, 0.0, 3.02]),
                      'y': numpy.array([-2.04, 0.0, 2.04])}

        #t & z are arbitrary and dx and dy seem reasonable
        #pylint: disable=invalid-name
        def image_sources(dx_dy_z_fluxbkp, t):
            """Return 3 overlapping sources with slightly different PSFs."""


            values = numpy.zeros((3, 3))
            d_dx = numpy.zeros((3, 3))
            d_dy = numpy.zeros((3, 3))
            d2_dxdy = numpy.zeros((3, 3))

            values[1, 1] = 10.0

            sources = []
            for dx, dy, z, flux_backup in dx_dy_z_fluxbkp:
                d_dx[1, 1] = ((1.0 - t) * dx / 150.0
                              +
                              (1.0 - t) * dy / 300.0)
                d_dy[1, 1] = (1.0 - t) * dx / 150.0 + 0.01 * z
#               print("Dx = " + repr(d_dx[1, 1]) + ", Dy = " + repr(d_dy[1, 1]))
                sources.append(
                    {
                        'x': 30.0 + dx,
                        'y': 30.0 + dy,
                        'z': z,
                        't': t,
                        'flux_backup': flux_backup,
                        'psf_args': {
                            'psf_parameters': {
                                'values': numpy.copy(values),
                                'd_dx': numpy.copy(d_dx),
                                'd_dy': numpy.copy(d_dy),
                                'd2_dxdy': numpy.copy(d2_dxdy),
                            },
                            'boundaries': boundaries
                        }
                    }
                )

            return sources
        #pylint: enable=invalid-name

        sources = [
            image_sources(
                [
                    (-17.3, -15.3, numpy.pi, None),
                    (12.3, -17.0, numpy.e, None),
                    (-15.0, 13.0, 1.0 / numpy.pi, None),
                    (17.3, 17.0, 1.0 / numpy.e, None)
                ],
                1.0
            ),
            image_sources(
                [
                    (0.0, -1.5, 1.0, None),
                    (-4.0, 1.5, 2.0, None),
                    (4.0, 1.5, 3.0, None)
                ],
                2.0
            ),
            image_sources(
                [
                    (5.5, -2.5, numpy.pi * numpy.e, None),
                    (5.5, 2.5, numpy.pi**2, None),
                    (0.0, 0.0, numpy.e**2, None),
                    (-5.5, -2.5, numpy.pi / numpy.e, None),
                    (-5.5, 2.5, -1.0, None)
                ],
                3.0
            ),
            image_sources(
                [
                    (-12.0, round(-5.0 * numpy.pi, 6), 10.0, None),
                    (-16.0, round(-5.0 * numpy.pi, 6), 10.0, None),
                    (14.5, -15.0, -10.0, None),
                    (16.2, -15.0, -10.0, None),
                    (-15.5, 13.15, 0.0, None),
                    (15.5, 14.83, -10.0, None),
                    (17.8, 15.0, 0.0, None),
                    (13.0, 15.8, -10.0, None),
                    (-17.8, -15.8, 0.0, None)
                ],
                4.0
            ),
            image_sources(#PROBLEM IMAGE
                [
                    (
                        dx,
                        dy,
                        1.0,
                        (
                            None if (
                                dx < -12.5
                                or
                                dx > 12.5
                                or
                                dy < -12.5
                                or
                                dy > 12.5
                            ) else (
                                10.0 * boundaries['x'][-1] * boundaries['y'][-1]
                            )
                        )
                    )
                    for dx in [-20.2, -15.3, -10.4, -5.5, 0.6, 5.7, 10.8,
                               15.9, 20.0]
                    for dy in [-20.2, -15.3, -10.4, -5.5, 0.6, 5.7, 10.8,
                               15.9, 20.0]
                ],
                5.0
            ),
            image_sources(
                [
                    (
                        dx,
                        dy,
                        3.0,
                        None
                    )
                    for dx in [-22.8, -15.2, -7.6, 0.0, 7.6, 15.2, 22.8]
                    for dy in [-22.8, -11.6, 0.0, 11.6, 22.8]
                ],
                5.0
            )
        ]
        self.run_test(sources=sources,
                      psffit_terms=['1', 'x', 'y', 't', 'x*t', 'y*t', 'z'])



if __name__ == '__main__':
    unittest.main(failfast=True)
