#!/usr/bin/env python3

"""Create noiseless FITS files to test various tools on."""

from subprocess import call
import unittest

import os
import os.path
import sys
import numpy
import h5py
from astropy.io import fits as pyfits

module_path = os.path.abspath(os.path.dirname(__file__))

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            module_path,
            '..',
            '..'
        )
    )
)

#Needs to be after os.path and sys to allow adding the seach path.
#pylint: disable=wrong-import-position

from astrowisp.fake_image.piecewise_bicubic_psf import PiecewiseBicubicPSF
from tests.fitpsf.utils import make_image_and_source_list
#pylint: enable=wrong-import-position


class TestPiecewiseBicubicNoiseless(unittest.TestCase):
    """Test piecewise bicubic PSF fitting on noiseless images."""

    fitpsf_executable = os.path.abspath(
        os.path.join(
            module_path,
            '..',
            '..',
            '..',
            'build',
            'exe',
            'fitpsf',
            'debug',
            'fitpsf'
        )
    )
#    fitpsf_executable = (
#        '/home/kpenev/projects/svn/HATpipe/source/'
#        'subpixel_sensitivity/src/build/exe/fitpsf/debug/fitpsf'
#    )

    #TODO: consider splitting into several functions
    #pylint: disable=too-many-locals
    def check_results(self, psf_fit_fname, sources, extra_variables):
        """
        Assert that fitted PSF map evaluates to expected PSFs for sources.

        Args:
            psf_fit_fname:    The name of the file produced by fitpsf to check
                the PSF map of.

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

        psf_fit_file = h5py.File(psf_fit_fname, 'r')
        map_group = psf_fit_file['PSFFit/Map']
        variables = {
            var: val[enabled_sources]
            for var, val in zip(['x', 'y'] + extra_variables,
                                map_group['Variables'][:])
        }

        psffit_terms = map_group.attrs['Terms'][0].decode()
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
        coefficients = map_group['Coefficients'][:]

        #Indices are: source index, variable, y boundary ind, x boundary ind
        fit_params = numpy.tensordot(term_list, coefficients, [1, 3])
        self.assertEqual(
            fit_params.shape,
            (num_sources, 4, num_x_boundaries, num_y_boundaries)
        )
        fluxes = psf_fit_file['PSFFit/Flux'][:]

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
    #pylint: enable=too-many-locals

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

        def grid_boundary_str(boundaries):
            """Return a comma separated list of the given grid boundaries."""

            return ','.join(str(b) for b in boundaries)

        if extra_variables is None:
            extra_variables = []
        for subpix_map in [numpy.ones((1, 1)),
                           numpy.ones((2, 2)),
                           numpy.array([[1.99, 0.01], [0.01, 1.99]]),
                           numpy.array([[0.5, 0.5], [0.5, 2.5]]),
                           numpy.array([[1.9], [0.1]]),
                           numpy.array([[2.0, 0.0], [0.0, 2.0]]),
                           numpy.array([[0.0, 0.0], [0.0, 4.0]])]:
            fname_start = os.path.join(module_path,
                                       'test_data',
                                       'noiseless_bicubic_psf.')
            filenames = dict(
                source_list=(fname_start + 'srclist'),
            )

            files_to_cleanup = [filenames['source_list']]

            if os.path.exists(files_to_cleanup[-1]):
                os.remove(files_to_cleanup[-1])

            for image_index, image_sources in enumerate(sources):
                filenames['image'] = (fname_start
                                      +
                                      str(image_index)
                                      +
                                      '.fits')
                filenames['psf_fit'] = (fname_start
                                        +
                                        str(image_index)
                                        +
                                        '.hdf5')
                print('Image sources:\n' + repr(image_sources))
                make_image_and_source_list(
                    sources=[dict(x=src['x'],
                                  y=src['y'],
                                  psf=PiecewiseBicubicPSF(**src['psf_args']),
                                  **{var: src[var] for var in extra_variables})
                             for src in image_sources],
                    extra_variables=extra_variables,
                    subpix_map=subpix_map,
                    filenames=filenames
                )
                if os.path.exists(filenames['psf_fit']):
                    os.remove(filenames['psf_fit'])

                files_to_cleanup.extend([filenames['image'],
                                         filenames['psf_fit']])


            subpix_fname = (fname_start + 'subpix.fits')
            pyfits.HDUList(
                [pyfits.PrimaryHDU(subpix_map)]
            ).writeto(
                subpix_fname,
                overwrite=True
            )
            files_to_cleanup.append(subpix_fname)

            with open(fname_start + 'cfg', 'w') as config:
                with open(
                    os.path.join(module_path, 'test_data', 'config_template.cfg'),
                    'r'
                ) as config_template:
                    config.write(
                        config_template.read()
                        %
                        dict(
                            source_list_fname=filenames['source_list'],
                            input_columns=','.join(['ID', 'x', 'y']
                                                   +
                                                   extra_variables),
                            terms=psffit_terms,
                            grid=(
                                grid_boundary_str(
                                    sources[0][0]['psf_args']['boundaries']['x']
                                )
                                +
                                ';'
                                +
                                grid_boundary_str(
                                    sources[0][0]['psf_args']['boundaries']['y']
                                )
                            ),
                            subpix_fname=subpix_fname
                        )
                    )
                files_to_cleanup.append(config.name)

            self.assertEqual(
                call([self.fitpsf_executable, '-c', files_to_cleanup[-1]]),
                0
            )

            for image_index, image_sources in enumerate(sources):
                self.check_results(
                    fname_start + str(image_index) + '.hdf5',
                    image_sources,
                    extra_variables
                )

            map(os.remove, files_to_cleanup)

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

    def test_isolated_sources(self):
        """Test fitting an image containing 8 well isolated sources."""

        psf_parameters = dict(values=numpy.zeros((3, 3)),
                              d_dx=numpy.zeros((3, 3)),
                              d_dy=numpy.zeros((3, 3)),
                              d2_dxdy=numpy.zeros((3, 3)))
        boundaries = dict(x=numpy.array([-2.0, 0.0, 2.0]),
                          y=numpy.array([-1.4, 0.0, 1.4]))

        sources = []

        psf_parameters['values'][1, 1] = 1.0
        sources.append(dict(x=15.0,
                            y=15.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 1.0
        sources.append(dict(x=45.0,
                            y=15.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(dict(x=15.0,
                            y=45.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.75
        psf_parameters['d_dy'][1, 1] = 1.00
        sources.append(dict(x=37.5,
                            y=45.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.5
        psf_parameters['d_dy'][1, 1] = 0.0
        sources.append(dict(x=30.0,
                            y=15.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.0
        psf_parameters['d_dy'][1, 1] = 0.5
        sources.append(dict(x=15.0,
                            y=30.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 0.5
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(dict(x=30.0,
                            y=45.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 0.5
        sources.append(dict(x=45.0,
                            y=30.0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))

        self.run_test(sources=[sources],
                      psffit_terms='{1, x, y}')

    def test_two_overlapping_sources(self):
        """Test fitting an image containing 2 sources all overlapping."""

        psf_args = dict(psf_parameters=dict(values=numpy.zeros((3, 3)),
                                            d_dx=numpy.zeros((3, 3)),
                                            d_dy=numpy.zeros((3, 3)),
                                            d2_dxdy=numpy.zeros((3, 3))),
                        boundaries=dict(x=numpy.array([-1.0, 0.0, 1.0]),
                                        y=numpy.array([-1.4, 0.0, 1.4])))

        psf_args['psf_parameters']['values'][1, 1] = 1.0

        sources = [dict(x=14.53,
                        y=14.02,
                        psf_args=psf_args),
                   dict(x=16.51,
                        y=14.03,
                        psf_args=psf_args)]
        self.run_test(sources=[sources], psffit_terms='{1}')

    def test_four_overlapping_sources(self):
        """Test fitting an image containing 4 overlapping sources."""

        psf_parameters = dict(
            values=numpy.zeros((3, 3)),
            d_dx=numpy.zeros((3, 3)),
            d_dy=numpy.zeros((3, 3)),
            d2_dxdy=numpy.zeros((3, 3))
        )

        boundaries = dict(
            x=numpy.array([-2.2, 0.0, 2.2]),
            y=numpy.array([-1.4, 0.0, 1.4])
        )

        psf_parameters['values'][1, 1] = 1.0

        sources = [dict(x=12.5,
                        y=13.5,
                        psf_args=dict(psf_parameters=dict(psf_parameters),
                                      boundaries=boundaries))]


        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        sources.append(dict(x=16.5,
                            y=13.5,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))

        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(dict(x=12.5,
                            y=15.5,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))

        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 1.0
        sources.append(dict(x=16.5,
                            y=15.5,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))

        self.run_test(sources=[sources], psffit_terms='{1, x*x, y*y}')

        for src in sources:
            src['enabled'] = 1
        psf_parameters['d_dx'] = numpy.zeros((3, 3))
        psf_parameters['d_dy'] = numpy.zeros((3, 3))
        psf_parameters['d2_dxdy'] = numpy.zeros((3, 3))
        psf_parameters['d_dx'][1, 1] = 1.0
        psf_parameters['d_dy'][1, 1] = 2.0
        psf_parameters['d2_dxdy'][1, 1] = 3.0
        sources.append(dict(x=4.0,
                            y=4.0,
                            enabled=0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        sources.append(dict(x=25.0,
                            y=25.0,
                            enabled=0,
                            psf_args=dict(psf_parameters=dict(psf_parameters),
                                          boundaries=boundaries)))
        self.run_test(sources=[sources],
                      psffit_terms='{1, x*x, y*y}',
                      extra_variables=['enabled'])

    def test_multi_image_with_extra_var(self):
        """Test fitting a series of 5 images and non-position variables."""

        boundaries = dict(x=numpy.array([-3.02, 0.0, 3.02]),
                          y=numpy.array([-2.04, 0.0, 2.04]))

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
#                print("Dx = " + repr(d_dx[1, 1]) + ", Dy = " + repr(d_dy[1, 1]))
                sources.append(
                    dict(
                        x=30.0 + dx,
                        y=30.0 + dy,
                        z=z,
                        t=t,
                        flux_backup=flux_backup,
                        psf_args=dict(
                            psf_parameters=dict(
                                values=numpy.copy(values),
                                d_dx=numpy.copy(d_dx),
                                d_dy=numpy.copy(d_dy),
                                d2_dxdy=numpy.copy(d2_dxdy),
                            ),
                            boundaries=boundaries
                        )
                    )
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
                    (4.0, 1.5, 3.0, None)],
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
            image_sources(
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
                      psffit_terms='{1, x, y, t, x*t, y*t, z}',
                      extra_variables=['t', 'z'])

if __name__ == '__main__':
    unittest.main(failfast=True)
