#!/usr/bin/env python3

"""Test SuperPhot's fit_star_shape module."""

import sys
import os.path
import unittest
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
from superphot import FitStarShape
from tests.utilities import FloatTestCase
#pylint: enable=wrong-import-position

class TestPSFFittingNoiseless(FloatTestCase):
    """Test piecewise bicubic PSF fitting on noiseless images."""

    def check_results(self, psf_map, sources, extra_variables):
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


