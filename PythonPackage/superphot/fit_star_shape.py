#!/usr/bin/env python3
"""Define the :class:`FitStarShape` class, which performs PSF/PRF fitting."""

from numbers import Number
from ctypes import\
    c_bool,\
    POINTER,\
    c_double,\
    c_char_p,\
    c_char,\
    c_ulong,\
    c_int
import numpy

from superphot._initialize_library import superphot_library
from superphot.superphot_io_tree import SuperPhotIOTree

class FitStarShape:
    """
    Fit for the PSF/PRF of stars and their flux.

    The PSF and PRF describe the distribution of light from a point source at
    infinity on the imaging device pixels.

        * PSF is the distribution of light produced by the optical system in the
          plane of the detector for a point source at infinity.

        * PRF is the response of a detector pixel at some offset from the
          center of the source. That is, the PRF is the PSF convolved with
          the sensitivity of detector pixels.

    Both representations use a general piece-wise polynomial model in which the
    area around the source location is split into a rectangular, but in general
    irregular, grid of cells. Over each cell the PSF/PRF is modeled as a
    bi-cubic polynomial. The resulting function is constraned to have continuous
    values, first order derivatives and x-y cross-derivative across cell
    boundaries. Further, it's value, first derivative and the x-y cross
    derivative are contsrainted to go to zero on the outside boundaries of the
    grid.

    .. _attributes:

    Attributes:
        _library_psf_fitter:    The library object for carrying out PSF/PRF
            fitting.

        _library_configuration:    Library configuration object set per the
            current :attr:`configuration`

        _result_tree:    The SuperPhotIOTree instance containing the last
            fittintg results, on None, if no fitting has been performed yet.

        mode(str):    Are we doing 'PSF' or 'PRF' fitting (case
            insensitive).

        configuration (dict):    The configuraiton for how to carry out PSF/PRF
            fitting. The following keys are used (others are ignored by this
            class):

            mode (str):
                What kind of fitting to do 'PSF' or 'PRF' (case insensitive).

            shape_terms (str):
                The terms the PSF is allowed to depend on. The EBNF grammar
                defining the language for this parameter is::

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

            grid (list of floats):
                A comma separated list of grid boundaries.  Can either be a
                single list, in which case it is used for both the horizontal
                and vertical boundaries. If different splitting is desired in
                the two directions, two lists should be supplied separated by
                ``;``. The first list should contain the vertical (x) boundaries
                and the second list gives the horizontal (y) ones.

            initial_aperture (float):
                This aperture is used to derive an initial guess for the
                amplitudes of sources when fitting for a piecewise bicubic PSF
                model by doing aperture photometry assuming a perfectly flat
                PSF.

            subpixmap (2D numpy array):
                The sub-pixel map, for PSF fitting only.

            smoothing (float):
                How much smoothing penalty to impose when fitting the PSF.
                ``None`` for no smoothing. Value can be both positive and
                negative and will always result in smoothing (less for negative
                values).

            max_chi2 (float):
                The value of the reduced chi squared above which sources are
                excluded from the fit. This can indicate non-point sources or
                sources for which the location is wrong among ohter things.

            pixel_rejection_threshold (float):
                A number defining individual pixels to exclude from the PSF fit.
                Pixels with fitting residuals (normalized by the standard
                deviation) bigger than this value are excluded. If zero, no
                pixels are rejected.

            max_abs_amplitude_change (float):
                The absolute root of sum squares tolerance of the source
                amplitude changes in order to declare the piecewise bicubic PSF
                fitting converged.

            max_rel_amplitude_change (float):
                The relative root of sum squares tolerance of the source
                amplitude changes in order to declare the piecewise bicubic PSF
                fitting converged.

            min_convergence_rate (float):
                If the rate of convergence falls below this threshold,
                iterations are stopped. The rate is calculated as the fractional
                decrease in the difference between the amplitude change and the
                value when it would stop, as determined by the
                :attr:`max_abs_amplitude_change` and
                :attr:`max_rel_amplitude_change` attributes.

            max_iterations (int):
                No more than this number if iterations will be performed. If
                convergence is not achieved before then, the latest estimates
                are output and an exception is thrown. A negative value allows
                infinite iterations. A value of zero, along with an initial
                guess for the PSF causes only the amplitudes to be fit for PSF
                fitting photometry with a known PSF. It is an error to pass a
                value of zero for this option and not specify and initial guess
                for the PSF.

            gain (float):
                The gain in electrons per ADU to assume for the input images.

            cover_grid (bool):
                If this option is true, all pixels that at least partially
                overlap with the grid are assigned to the corresponding source.
                This option is ignored for sdk PSF models.

            src_min_signal_to_noise (float):
                How far above the background (in units of RMS) should pixels be
                to still be considered part of a source. Ignored if the
                piecewise bibucic PSF grid is used to select source pixels
                (cover-bicubic-grid option).

            src_max_aperture (float):
                If this option has a positive value, pixels are assigned to
                sources in circular apertures (the smallest such that all pixels
                that pass the signal to noise cut are still assigned to the
                source). If an aperture larger than this value is required, an
                exception is thrown.

            src_max_sat_frac (float):
                If more than this fraction of the pixels assigned to a source
                are saturated, the source is excluded from the fit.

            src_min_pix (int):
                The minimum number of pixels that must be assigned to a source
                in order to include the source is the PSF fit.

            src_max_pix (int):
                The maximum number of pixels that car be assigned to a source
                before excluding the source from the PSF fit.

            src_max_count (int):
                The maximum number of sources to include in the fit for the PSF
                shape. The rest of the sources get their amplitudes fit and are
                used to determine the overlaps. Sources are ranked according to
                the sum of (background excess)^2/(pixel variance+background
                variance) of their individual non-saturated pixels.

            bg_min_pix (int):
                The minimum number of pixels a background estimate must be based
                on in order to include the source in shape fitting.

            magnitude_1adu (float):
                The magnitude that corresponds to a flux of 1ADU.


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

    _default_configuration = dict(subpixmap=numpy.ones((1, 1), dtype=c_double),
                                  smoothing=None,
                                  max_chi2=100.0,
                                  pixel_rejection_threshold=100.0,
                                  max_abs_amplitude_change=0.0,
                                  max_rel_amplitude_change=1e-6,
                                  min_convergence_rate=-numpy.inf,
                                  max_iterations=1000,
                                  gain=1.0,
                                  cover_grid=True,
                                  src_min_signal_to_noise=3.0,
                                  src_max_aperture=10.0,
                                  src_max_sat_frac=1.0,
                                  src_min_pix=5,
                                  src_max_pix=1000,
                                  src_max_count=10000,
                                  bg_min_pix=50,
                                  magnitude_1adu=10.0)

    #Many return statements make sense in this case.
    #pylint: disable=too-many-return-statements
    @staticmethod
    def _format_config(param_value):
        """Format config param for passing to SuperPhot PSF fitting lib."""

        prefix = b''
        if param_value[0].startswith('src_'):
            param_value = (param_value[0][4:], param_value[1])
            prefix = b'src.'
        if param_value[0].startswith('bg_'):
            param_value = (param_value[0][3:], param_value[1])
            prefix = b'bg.'

        elif param_value[0] == 'cover_grid':
            return (b'src.cover-bicubic-grid',
                    repr(param_value[1]).encode('ascii'))
        elif param_value[0] == 'grid':
            grid = param_value[1]
            return (
                b'psf.bicubic.grid',
                (
                    ','.join(map(str, grid)) if isinstance(grid[0], Number)
                    else ';'.join([','.join(map(repr, grid_part))
                                   for grid_part in grid])
                ).encode('ascii')
            )
        elif param_value[0] == 'subpixmap':
            return ()
        elif param_value[0] == 'smoothing' and param_value[1] is None:
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
    #pylint: enable=too-many-return-statements

    def __init__(self,
                 *,
                 mode,
                 shape_terms,
                 grid,
                 initial_aperture,
                 **other_configuration):
        """
        Set-up an object ready to perform PSF/PRF fitting.

        Args:
            All the configuration attributes of the class can be configured by
            passing them as keyword arguments.

        Returns:
            None
        """

        self.mode = mode.upper()
        assert self.mode in ['PSF', 'PRF']
        self.configuration = dict(self._default_configuration)
        self.configuration.update(shape_terms=shape_terms,
                                  grid=grid,
                                  initial_aperture=initial_aperture,
                                  **other_configuration)

        self._library_configuration = (
            superphot_library.create_psffit_configuration()
        )
        self.configure()

        self._result_tree = None

    def configure(self, **configuration):
        """
        Modify the currently defined configuration.

        Args:
            **configuration:    See :attr:`configuration`\ .

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
            map(self._format_config, self.configuration.items()),
            (
                c_bool(self.mode == 'PRF'),
                self._library_configuration,
                b'psf.model',
                b'bicubic'
            )
        ) + (b'',)
        superphot_library.update_psffit_configuration(*config_arguments)

    def fit(self, image_sources, backgrounds):
        """
        Fit for the shape of the sources in a collection of imeges.

        Args:
            image_sources (list of 4-tuples):    Each entry consists of:

                1. The pixel values of the calibratred image

                2. The error estimates of the pixel values

                3. Mask flags of the pixel values.

                4. Sources to process, defining at least the following
                   quantities:

                       * **ID** (string): some unique identifier for the source

                       * **x** (float): The x coordinate of the source
                         center in pixels

                       * **y** (float): See ``x``

                   May define additional quantities on which the PSF shape is
                   allowed to depend.
                   
                   The source list can be either a numy record array with field
                   names as keys or a dictionary with field names as keys and
                   1-D numpy arrays of identical lengths as values.

            backgrounds:    The measured backgrounds under the sources (instance
                of :class:`BackgroundExtractor`\ ).

        Returns:
            None. Use :meth:`get_last_fit_result` to obtain the results.
        """

        def create_image_arguments():
            """
            Create the three image arguments for piecewise_bicubic_fit.

            Args:
                None

            Returns:
                tuple:
                    POINTER(POINTER(c_double)):
                        The pixel_values argument to the piecewise_bicubic_fit
                        library function

                    POINTER(POINTER(c_double)):
                        The pixel_errors argument to the piecewise_bicubic_fit
                        library function

                    POINTER(POINTER(c_char)):
                        The pixel_masks argument to the piecewise_bicubic_fit
                        library function

                    int:
                        The number of images to simultaneously process.

                    int:
                        The common x resolution of the images.

                    int:
                        The common y resolution of the images.

                Raises:
                    AssertionError:    If the shapes of the images do not all
                        match.
            """

            number_images = len(image_sources)
            image_y_resolution, image_x_resolution = image_sources[0][0].shape

            for entry in image_sources:
                for image in entry[:3]:
                    assert image.shape == (image_y_resolution,
                                           image_x_resolution)

            return (
                (POINTER(c_double) * number_images)(
                    *(
                        entry[0].ctypes.data_as(POINTER(c_double))
                        for entry in image_sources
                    )
                ),
                (POINTER(c_double) * number_images)(
                    *(
                        entry[1].ctypes.data_as(POINTER(c_double))
                        for entry in image_sources
                    )
                ),
                (POINTER(c_char) * number_images)(
                    *(
                        entry[2].ctypes.data_as(POINTER(c_char))
                        for entry in image_sources
                    )
                ),
                number_images,
                image_x_resolution,
                image_y_resolution
            )

        def get_column_names():
            """
            Return the list of columns defined for the input sources.

            Args:
                None

            Returns:
                list:
                    The column names which participate in the PSF expansion.

            Raises:
                AssertionError:    If the columns defined for all images do not
                    match or if any of the minimum required columns: 'ID', 'x',
                    'y' is missing.
            """

            if isinstance(image_sources[0][3], numpy.ndarray):
                column_names = image_sources[0][3].dtype.names
                for entry in image_sources:
                    assert entry[3].dtype.names == column_names
            else:
                column_names = image_sources[0][3].keys()
                column_name_set = set(column_names)
                for entry in image_sources:
                    assert set(entry[3].keys()) == column_name_set

            column_names = list(column_names)

            assert 'ID' in column_names
            assert 'x' in column_names
            assert 'y' in column_names


            column_names.remove('ID')

            return column_names

        def create_column_data(column_names):
            """
            Create the column_data array required by create_source_arguments.

            Args:
                column_names([str]):    The columns, other than ID defined for
                    the input list of sources.

            Returns:
                [numpy.ndarray]:
                    See column_data argument of create_source_arguments.
            """

            column_data = [
                numpy.empty((len(column_names), len(entry[3]['ID'])),
                            dtype=c_double)
                for entry in image_sources
            ]
            for image_index, entry in enumerate(image_sources):
                for column_index, column_name in enumerate(column_names):
                    column_data[image_index][column_index, :] = (
                        entry[3][column_name]
                    )
            return column_data

        def create_source_arguments(column_names, column_data):
            """
            Create the arguments defining the sources for piecewise_bicubic_fit.

            Args:
                column_names([str]):    The columns, other than ID defined for
                    the input list of sources.

                column_data([numpy.ndarray]):    List of one 2-D array for each
                    image with the first array index going accross columns and
                    the second array index going accross sources.

            Returns:
                tuple:
                    POINTER(c_char_p):    The column_names argument to the
                        piecewise_bicubic_fit library function.

                    POINTER(POINTER(c_char_p)):    The source_ids argument to
                        the piecewise_bicubic_fit library function.

                    POINTER(POINTER(c_double)):    The column_data argument to
                        the piecewise_bicubic_fit library function.

                    numpy.array(c_ulong):    1-D array contining the number of
                        sources in each image.

                    int:    The number of columns in the column_data.
            """

            number_images = len(image_sources)
            number_columns = len(column_names)

            return (
                (c_char_p * number_columns)(
                    *(
                        c_char_p(colname.encode('ascii'))
                        for colname in column_names
                    )
                ),
                (POINTER(c_char_p) * number_images)(
                    *(
                        (c_char_p * len(entry[3]['ID']))(
                            *(
                                (
                                    source_id if isinstance(source_id, bytes)
                                    else source_id.encode('ascii')
                                )
                                for source_id in entry[3]['ID']
                            )
                        )
                        for entry in image_sources
                    )
                ),
                (POINTER(c_double) * number_images)(
                    *(
                        columns.ctypes.data_as(POINTER(c_double))
                        for columns in column_data
                    )
                ),
                numpy.array([len(entry[3]['ID']) for entry in image_sources],
                            dtype=c_ulong),
                number_columns
            )

        column_names = get_column_names()
        column_data = create_column_data(column_names)
        result_tree = SuperPhotIOTree(self._library_configuration)
        if not superphot_library.piecewise_bicubic_fit(
                *create_image_arguments(),
                *create_source_arguments(column_names, column_data),
                (
                    len(backgrounds)
                    *
                    superphot_library.create_background_extractor.restype
                )(
                    *(bg.library_extractor for bg in backgrounds)
                ),
                self._library_configuration,
                self.configuration['subpixmap'],
                self.configuration['subpixmap'].shape[1],
                self.configuration['subpixmap'].shape[0],
                result_tree.library_tree
        ):
            raise RuntimeError("Star shape fitting failed to converge!")
        return result_tree

    def __del__(self):
        r"""Destroy the configuration object created in :meth:`__init__`\ ."""

        superphot_library.destroy_psffit_configuration(
            self._library_configuration
        )

if __name__ == '__main__':
    fitprf = FitStarShape(mode='prf',
                          shape_terms='{1}',
                          grid=[-1.0, 0.0, 1.0],
                          initial_aperture=2.0,
                          smoothing=None,
                          min_convergence_rate=0.0)

    #Debugging code
    #pylint: disable=protected-access
    tree = SuperPhotIOTree(fitprf._library_configuration)
    #pylint: enable=protected-access
    print('BG tool: ' + repr(tree.get('bg.tool', str)))
    print('PSF terms: ' + repr(tree.get('psffit.terms', str)))
    print('Max chi squared: '
          +
          repr(tree.get('psffit.max_chi2', c_double)))
    print('Maximum iterations: '
          +
          repr(tree.get('psffit.max_iterations', c_int)))
    print('Pixel rejection threshold: '
          +
          repr(tree.get('psffit.pixrej', c_double)))
