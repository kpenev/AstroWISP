#!/usr/bin/env python3
"""Define the :class:`FitStarShape` class, which performs PSF/PRF fitting."""

from numbers import Number
from ctypes import cdll, c_void_p, c_bool
from ctypes.util import find_library
import numpy

from _initialize_core_library import _initialize_core_library

#Naming convention imitates the one by ctypes.
#pylint: disable=invalid-name
#Type checking place holders require no content.
#pylint: disable=too-few-public-methods

class _c_fitting_configuration(c_void_p):
    """Placeholder for the FittingConfiguration opaque struct of fitpsf lib."""

#pylint: enable=invalid-name
#pylint: enable=too-few-public-methods

def _initialize_library():
    """Prepare the SuperPhot fitpsf library for use."""

    fitpsf_library_fname = find_library('superphotfitpsf')
    psf_library_fname = find_library('superphotpsf')
    io_library_fname = find_library('superphotio')
    if fitpsf_library_fname is None:
        raise OSError('Unable to find the SuperPhot PSF fitting library.')
    if psf_library_fname is None:
        raise OSError('Unable to find the SuperPhot PSF library.')
    if io_library_fname is None:
        raise OSError('Unable to find the SuperPhot I/O library.')
    cdll.LoadLibrary(io_library_fname)
    cdll.LoadLibrary(psf_library_fname)
    library = cdll.LoadLibrary(fitpsf_library_fname)

    _initialize_core_library(library)

    library.create_psffit_configuration.argtypes = []
    library.create_psffit_configuration.restype = _c_fitting_configuration

    library.destroy_psffit_configuration.argtype = [
        library.create_psffit_configuration.restype
    ]

    library.update_psffit_configuration.restype = None

    return library

class FitStarShape:
    """
    Fit for the PSF/PRF of stars and their flux.

    .. _attributes:

    Attributes:
        _library_psf_fitter:    The library object for carrying out PSF/PRF
            fitting.

        _library_configuration:    Library configuration object set per the
            current :attr:`configuration`

        _library_subpixmap:    Library sub-pixel sensitivity map object matching
            the current configuration.

        mode(str):    Are we doing `'PSF'` or `'PRF'` fitting (case
            insensitive).

        configuration (dict):    The configuraiton for how to carry out PSF/PRF
            fitting. See the keyword arguments to :meth:`__init__`.

    Example:
        Create and configure a PRF fitting object, allowing up to third order
        dependence on image position, on a grid which splits the area around the
        source in 16 squares of 2pix by 2pix size each and using an aperture
        with 5 pixel radius for the initial estimate of source amplitudes:

        >>> from superphot import FitStarShape
        >>> fitprf = FitStarShape(mode='prf',
        >>>                       shape_terms='O3{x, y}',
        >>>                       grid=[-4.0, -2.0, 0.0, 2.0, 4.0],
        >>>                       initial_aperture=5.0)
    """

    _library = _initialize_library()

    @staticmethod
    def _format_config(param_value):
        """Format config param for passing to SuperPhot PSF fitting lib."""

        prefix = b''
        if param_value[0] == 'grid':
            grid = param_value[1]
            return (
                b'grid',
                (
                    ','.join(map(str, grid)) if isinstance(grid[0], Number)
                    else ';'.join([','.join(map(repr, grid_part))
                                   for grid_part in grid])
                ).encode('ascii')
            )
        elif param_value[0] == 'subpixmap':
            return ()
        elif param_value[0] == 'shape_terms':
            return b'psf.terms', param_value[1].encode('ascii')
        elif param_value[0] == 'pixel_rejection_threshold':
            return (b'psf.bicubic.pixrej',
                    repr(param_value[1]).encode('ascii'))
        elif param_value[0] in ['max_iterations',
                                'max_chi2',
                                'min_convergence_rate']:
            prefix = b'psf.'
        elif param_value[0] in ['max_abs_amplitude_change',
                                'max_rel_amplitude_change',
                                'initial_aperture',
                                'grid',
                                'smoothing']:
            prefix = b'psf.bicubic.'
        return (
            prefix + param_value[0].replace('_', '-').encode('ascii'),
            (
                param_value[1] if isinstance(param_value[1], str)
                else repr(param_value[1])
            ).encode('ascii')
        )

    def __init__(self,
                 *,
                 mode,
                 shape_terms,
                 grid,
                 initial_aperture,
                 subpixmap=numpy.ones((1, 1), dtype=float),
                 smoothing=None,
                 max_chi2=100.0,
                 pixel_rejection_threshold=100.0,
                 max_abs_amplitude_change=0.0,
                 max_rel_amplitude_change=1e-6,
                 min_convergence_rate=-numpy.inf,
                 max_iterations=1000):
        """
        Set-up an object ready to perform PSF/PRF fitting.

        Args:
            mode (str):    What kind of fitting to do 'PSF' or 'PRF' (case
                insensitive).

            shape_terms (str):    The terms the PSF is allowed to depend on. The
                EBNF grammar defining the language for this parameter is::

                    (* items in angle brackets (< or >) are assumed to be     *)
                    (* obvious and thus are not defined.                      *)

                    termchar = <ascii character> - "," - "}" ;

                    (* mathematical expressions  involving variables,         *)
                    (* floating point numbers and pi. The complete list of    *)
                    (* mathematical functions from c++99's cmath library are  *)
                    (* supported.                                             *)
                    term = termchar , { termchar } ;


                    (* simple listing of terms to include.                    *)
                    list = "{" , term , { "," , term } , "}" ;

                    (* Expands to all polynomial terms of up to combined      *)
                    (* order <integer> of the entries in list                 *)
                    poly = "O" , <integer> , list ;

                    set = list | poly ;


                    (* expands to the cross product of all sets.              *)
                    cross = set , { "*" , set } ;


                    (* merge the terms of all cross products together.        *)
                    expression = cross , { \"+\" , cross } ;

            grid (list of floats):    A comma separated list of grid boundaries.
                Can either be a single list, in which case it is used for both
                the horizontal and vertical boundaries. If different splitting
                is desired in the two directions, two lists should be supplied
                separated by ``;``. The first list should contain the vertical
                (x) boundaries and the second list gives the horizontal (y)
                ones.

            initial_aperture (float):    This aperture is used to derive an
                initial guess for the amplitudes of sources when fitting for a
                piecewise bicubic PSF model by doing aperture photometry
                assuming a perfectly flat PSF.

            subpixmap (2D numpy array):    The sub-pixel map, for PSF fitting
                only.

            smoothing (float):    How much smoothing penalty to impose when
                fitting the PSF.  ``None`` for no smoothing. Value can be both
                positive and negative and will always result in smoothing (less
                for negative values).

            max_chi2 (float):    The value of the reduced chi squared above
                which sources are excluded from the fit. This can indicate
                non-point sources or sources for which the location is wrong
                among ohter things.

            pixel_rejection_threshold (float):    A number defining individual
                pixels to exclude from the PSF fit.  Pixels with fitting
                residuals (normalized by the standard deviation) bigger than
                this value are excluded. If zero, no pixels are rejected.

            max_abs_amplitude_change (float):    The absolute root of sum
                squares tolerance of the source amplitude changes in order to
                declare the piecewise bicubic PSF fitting converged.

            max_rel_amplitude_change (float):    The relative root of sum
                squares tolerance of the source amplitude changes in order to
                declare the piecewise bicubic PSF fitting converged.

            min_convergence_rate (float):    If the rate of convergence falls
                below this threshold, iterations are stopped. The rate is
                calculated as the fractional decrease in the difference between
                the amplitude change and the value when it would stop, as
                determined by the :attr:`max_abs_amplitude_change` and
                :attr:`max_rel_amplitude_change` attributes.

            max_iterations (int):    No more than this number if iterations will
                be performed. If convergence is not achieved before then, the
                latest estimates are output and an exception is thrown. A
                negative value allows infinite iterations. A value of zero,
                along with an initial guess for the PSF causes only the
                amplitudes to be fit for PSF fitting photometry with a known
                PSF. It is an error to pass a value of zero for this option and
                not specify and initial guess for the PSF.

        Returns:
            None
        """

        self.mode = mode.upper()
        assert self.mode in ['PSF', 'PRF']
        self.configuration = dict(
            subpixmap=subpixmap,
            shape_terms=shape_terms,
            max_chi2=max_chi2,
            grid=grid,
            pixel_rejection_threshold=pixel_rejection_threshold,
            initial_aperture=initial_aperture,
            max_abs_amplitude_change=max_abs_amplitude_change,
            max_rel_amplitude_change=max_rel_amplitude_change,
            smoothing=smoothing,
            min_convergence_rate=min_convergence_rate,
            max_iterations=max_iterations
        )

        self._library_configuration = (
            self._library.create_psffit_configuration()
        )
        config_arguments = sum(
            map(self._format_config, self.configuration.items()),
            (c_bool(self.mode == 'PRF'), b'psf.model', b'bicubic')
        ) + (b'',)
        self._library.update_psffit_configuration(
            self._library_configuration,
            *config_arguments
        )

        self._library_subpixmap = None
        self.set_subpix_map(subpixmap)

    def configure(self, **configuration):
        """
        Modify the currently defined configuration.

        Args:
            **configuration:    See the keyword arguments of :meth:`__init__`.

        Returns:
            None
        """

        for k in configuration:
            if k not in self.configuration:
                raise KeyError('Unrecognized configuration parameter: '
                               +
                               repr(k))

        if 'mode' in configuration:
            self.mode = configuration['mode'].upper()
            assert self.mode in ['PSF', 'PRF']

        self.configuration.update(configuration)

        config_arguments = sum(
            map(self._format_config, configuration.items()),
            (c_bool(self.mode == 'PRF'),)
        ) + (b'',)
        self._library.update_psffit_configuration(
            self._library_configuration,
            *config_arguments
        )

        if 'subpixmap' in configuration:
            self.set_subpix_map(configuration['subpixmap'])

    def set_subpix_map(self, subpixmap):
        """
        Modify the sub-pixel sensitivity map to use for PSF fitting.

        Args:
            subpixmap:    See same name keyword argument to :meth:`__init__`

        Returns:
            None
        """

        self.configuration['subpixmap'] = subpixmap

        if self._library_subpixmap is not None:
            self._library.destroy_core_subpixel_map(self._library_subpixmap)

        self._library_subpixmap = self._library.create_core_subpixel_map(
            subpixmap.shape[1],
            subpixmap.shape[0],
            subpixmap
        )

    def __call__(self, image_sources):
        """
        Fit for the shape of the sources in a collection of imeges.

        Args:
            image_sources (dict):    Dictionary, with keys - the filename of the
                reduced frame (FITS) to fit the PSF/PRF of and values - list of
                sources to process, defining at least the following quantities:

                    * **ID** (string): some unique identifier for the source

                    * **x** (float): The x coordinate of the source center in pixels

                    * **y** (float): See ``x``

                May define additional quantities on which the PSF shape is
                allowed to depend. This can be either a numy record array with
                field names as keys or a dictionary with field names as keys
                and 1-D numpy arrays of identical lengths as values.

        Returns:
            dict:

                * keys:
                    The filenames of the reduced frames used in fitting.

                * values:
                    2-tuples containing the following arrays:

                    * 4D numpy array:
                        The coefficients of the PSF/PRF map.

                    * 1D numpy array:
                        The best-fit fluxes of the sources in the same order as
                        specified in the :obj:`image_sources` argument.
        """

        raise Exception('Not implemented')

    def __del__(self):
        r"""Destroy the configuration object created in :meth:`__init__`\ ."""

        self._library.destroy_psffit_configuration(self._library_configuration)

if __name__ == '__main__':
    fitprf = FitStarShape(mode='prf',
                          shape_terms='O3{x, y}',
                          grid=[-4.0, -2.0, 0.0, 2.0, 4.0],
                          initial_aperture=5.0,
                          smoothing=0.0,
                          min_convergence_rate=0.0)
