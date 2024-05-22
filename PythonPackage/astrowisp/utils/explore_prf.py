#!/usr/bin/env python3

#pylint: disable=too-many-lines

"""Create some plots helpful for picking PSF grid."""

import os.path
import subprocess
from ctypes import c_double, c_char

from matplotlib import pyplot
import scipy
import scipy.spatial
import scipy.stats
#import numpy
from scipy.interpolate import SmoothBivariateSpline
try:
    #pylint: disable=import-error
    import xalglib
    #pylint: enable=import-error
finally:
    pass
from astropy.io import fits
from configargparse import\
    ArgumentParser,\
    Action as ArgparseAction,\
    DefaultsFormatter


from general_purpose_python_modules.kelly_colors import kelly_colors

from astrowisp import BackgroundExtractor, FitStarShape, SubPixPhot
from astrowisp.utils.file_utilities import\
    get_fname_pattern_substitutions,\
    prepare_file_output
from astrowisp.utils import flux_from_magnitude

#pylint: disable=too-many-statements
def parse_command_line(parser=None,
                       assume_sources=False,
                       add_config_file=True,
                       add_frame_arg=True):
    """
    Return the command line arguments as attributes of an object.

    Args:
        parser(ArgumentParser):    If not None, should be a valid command line
            parser to which further arguments are addded by this method.

    Returns:
        See ArgumentParser.parse_args().
    """

    #Interface mandated by argparse.
    #pylint: disable=too-few-public-methods
    class ValidateBinning(ArgparseAction):
        """Properly parse the different type arguments specifying binning."""

        def __call__(self, parser, namespace, values, option_string=None):
            """Add dictionary with binning configuration to namespace."""

            assert len(values) == 2

            setattr(namespace,
                    self.dest,
                    {'statistic': values[0], 'bins': int(values[1])})
    #pylint: enable=too-few-public-methods

    def parse_image_split(split_str):
        """Parse the image split argument to a tuple of (direction, value)."""

        direction, value = split_str.split('=')
        direction = direction.strip()
        assert direction in ['x', 'y']

        return direction, int(value.strip())

    def parse_slice(slice_str):
        """Parse slice arguments to dictionaries directly usable as kwargs."""

        direction, slice_str = slice_str.split('=')
        direction = direction.strip()
        assert direction in ['x', 'y']

        offset, thickness = (float(val_str.strip())
                             for val_str in slice_str.split('+-'))
        return {direction + '_offset': offset,
                'thickness': thickness}

    if parser is None:
        parser = ArgumentParser(
            description=__doc__,
            default_config_files=['explore_prf.cfg'],
            formatter_class=DefaultsFormatter,
            ignore_unknown_config_file_keys=True
        )

    if add_frame_arg:
        parser.add_argument(
            'frame_fname',
            help='The full path of the FITS file to use for creating the plots.'
        )
    parser.add_argument(
        '--prf-range', '-r',
        default=(8.0, 8.0, 4.0, 4.0),
        nargs=4,
        type=float,
        help="Width, height, x and y offset of the source extraction center "
        "from the lower left corner. Default: `%(default)s'."
    )
    parser.add_argument(
        '--background-annulus', '--bg', '-b',
        type=float,
        nargs=2,
        default=(6.0, 7.0),
        help='The inner and outer radius of the annulus used to measure the '
        'background of sources. Default: %(default)s.'
    )
    parser.add_argument(
        '--flux-aperture', '--aperture', '-a',
        type=float,
        default=4.0,
        help='The aperture to use for measuring the flux (normalization of the '
        'pixel values when plotting. Default: %(default)s.'
    )
    parser.add_argument(
        '--slice', '-s',
        type=parse_slice,
        action='append',
        default=None,
        help='Add more slices to show. Each slice is specified as an offset '
        'along one if the demensions (x or y) and a range around that to '
        'include in the plot. White space around tokens is allowed and '
        'ignored. Example: "x = 0 +- 0.1". By default plot slices at x=0 and '
        'y=0 with width of 0.2.'
    )
    parser.add_argument(
        '--split-image',
        type=parse_image_split,
        action='append',
        default=[],
        help='Specify another boundary (in x or y) to spling the image into '
        'regions. A separate plot of the PRF is generated for each region. The '
        'format is like "x/y=<value>", and the option can be specified multiple'
        ' times.'
    )
    parser.add_argument(
        '--discard-image-boundary',
        action='store_true',
        help='If passed, slices on the outside boundary of the image are not '
        'plotted.'
    )
    parser.add_argument(
        '--error-scale',
        type=float,
        default=0.1,
        help='Scale the error bars by this value when plotting to make a more '
        'readable plot. Default: %(default)s.'
    )
    parser.add_argument(
        '--error-threshold',
        type=float,
        default=0.1,
        help='Points with error bars larger than this are not included in the '
        'plot. Intended to avoid introducing points that are too noisy. '
        'Default: %(default)s'
    )
    parser.add_argument(
        '--add-binned',
        nargs=2,
        default={},
        action=ValidateBinning,
        metavar=('STATISTIC', 'NBINS'),
        help='If supplied, in addition to pixel values, a binned curve is also '
        'displayed using the specified binning statistic and number of bins.'
    )
    if not assume_sources:
        parser.add_argument(
            '--trans-pattern', '-t',
            default=os.path.join('%(FITS_DIR)s',
                                 '..',
                                 'ASTROM',
                                 '%(FITS_ROOT)s.trans'),
            help="A pattern with substitutions involving any FITS header "
            "keyword, `'%%(FITS_DIR)s'` (directory containing the frame), "
            "and/or `'%%(FITS_ROOT)s'` (base filename of the frame without the "
            "`fits` or `fits.fz` extension) that expands to the `.trans` file "
            "corresponding to the input frame. Default: '%(default)s'."
        )
        parser.add_argument(
            '--catalogue-pattern',
            default=os.path.join('%(FITS_DIR)s',
                                 '..',
                                 'MASTERS',
                                 '%(FILTER)c_catalogue.ucac4'),
            help="A pattern with substitutions involving any FITS header "
            "keywordd, `'%%(FITS_DIR)s'` (directory containing the frame), "
            "and/or `'%%(FITS_ROOT)s'` (base filename of the frame without the "
            "`fits` or `fits.fz` extension) that expands to the full path of "
            "the catalogue containing all sources in the image. Default: "
            "'%(default)s'."
        )

    spline_config = parser.add_argument_group(
        title='Spline configuration',
        description='Options controlling the bicubic spline interpolation used '
        'to represent the PRF in each image piece.'
    )
    spline_config.add_argument(
        '--spline-method',
        choices=['alglib', 'scipy', 'none'],
        default='alglib',
        help='Which of the supported spline fitting methods to use. Default: '
        '%(default)s'
    )
    spline_config.add_argument(
        '--spline-resolution',
        nargs=2,
        type=int,
        default=(20, 20),
        help="The grid to use if the spline method is 'alglib'. Ignored if the "
        "method is 'scipy'. Default: %(default)s."
    )
    spline_config.add_argument(
        '--spline-smoothing',
        type=float,
        default=1.0,
        help='For scipy splines, the bi-cubic spline derived has a smoothing '
        'factor equal to this value times the number of points at which PRF '
        'measurements are available. For alglib splines, this is directly '
        'the smoothing penalty used. Default: %(default)s.'
    )
    spline_config.add_argument(
        '--spline-pad-fraction',
        type=float,
        default=0.01,
        help='The fraction of the PRF grid box to add with zero-valued padding '
        'past the edges. Default: %(default)s.'
    )
    spline_config.add_argument(
        '--spline-pad-npoints',
        type=int,
        default=100,
        help='The number of zero valued points, in each direction to add within'
        'each padding region. Default: %(default)s.'
    )
    spline_config.add_argument(
        '--spline-spacing',
        type=int,
        default=300,
        help='The number of evenly spaced points to pass to the spline for x '
        'and y. Default: %(default)s.'
    )

    plot_config = parser.add_argument_group(
        title='Plotting configuration',
        description='Options that control the layout and other aspects of the '
        'plots generated.'
    )
    plot_config.add_argument(
        '--figure-dpi',
        type=int,
        default=300,
        help='The resolution to set for generated plots.'
    )
    plot_config.add_argument(
        '--marker-size',
        type=float,
        default=2.0,
        help='The size of the markers to use in plots.'
    )
    plot_config.add_argument(
        '--plot-y-range',
        default=None,
        type=float,
        nargs=2,
        help='If specified, the generated plot displays exactly the given range'
        ' of y values.'
    )
    plot_config.add_argument(
        '--plot-3d-spline',
        action='store_true',
        help='If true, the plots will be generated as 3D images of the spline'
    )
    plot_config.add_argument(
        '--save-3d-spline-plot',
        default=None,
        help="This should be a template for the filename for plotting the "
        "3D prf involving `'%%(XSPLIT)s'` and `'%%(YSPLIT)s'` which get "
        "substituted by the image slices given, and/or any FITS header keyword,"
        " and/or `'%%(FITS_DIR)s'` (directory containing the frame), "
        "`'%%(FITS_ROOT)s'` (base filename of the frame without any "
        "extensions). This should be enabled with plot-3d-spline if you want "
        "plots to be saved otherwise plots will just be shown"
    )
    plot_config.add_argument(
        '--plot-3d-spline-contour-stride',
        default=60,
        type=int,
        help='The stride amount (int) for the 3D contour plot generated by '
        'plot-3d-spline'
    )
    plot_config.add_argument(
        '--plot-3d-spline-contour-linewidths',
        default=60.0,
        type=float,
        help='The The line width of the contour lines (float) for the 3D '
        'contour plot generated by plot-3d-spline'
    )
    plot_config.add_argument(
        '--plot-multi-image',
        action='store_true',
        help='If true, the plots will be generated as multiple separate images,'
             ' instead of stacked on top of one another'
    )
    plot_config.add_argument(
        '--plot-entire-prf',
        action='store_true',
        help='If true, will generate a plot of the entire prf on a 3d scale.'
    )
    plot_config.add_argument(
        '--save-prf-plot',
        default=None,
        help="This should be a template for the filename for plotting the "
        "entire prf involving `'%%(XSPLIT)s'` and `'%%(YSPLIT)s'` which get "
        "substituted by the image slices given, and/or any FITS header keyword,"
        " and/or `'%%(FITS_DIR)s'` (directory containing the frame), "
        "`'%%(FITS_ROOT)s'` (base filename of the frame without any "
        "extensions)."
    )
    plot_config.add_argument(
        '--save-plot-pattern',
        default=None,
        help="If this option is not specified, plots will be displayed, "
        "otherwise, this should be a template for the filename involving "
        "any FITS header keywords, `'%%(FITS_DIR)s'` (directory containing the "
        "frame), `'%%(FITS_ROOT)s'` (base filename of the frame without any "
        "extensions), %%(SLICEDIR), and %%(SLICEOFFSET) substitutions which get"
        " substituted by the direction (`\'x\'` or `\'y\'`) and the offset of "
        "the slice in the plot. If using plot-multi-image option include "
        "%%(XSPLIT)s and %%(YSPLIT)s which get substituted by the image "
        "slices given."
    )
    parser.add_argument(
        '--skip-existing-plots',
        action='store_true',
        help='If this argument is enabled, if the output file for a plot '
        'exists, the plot is not re-generated. This allows saving the time to '
        'calculate the required data if all plots exist.'
    )
    parser.add_argument(
        '--overwrite-plots',
        action='store_true',
        help='If a plot exists, without this argument or --skip-existing-plots,'
        ' the script exists with an error.'
    )

    if add_config_file:
        parser.add_argument(
            '--config', '-c',
            is_config_file=True,
            help='Config file to use instead of default.'
        )

    #TODO write common line argument for markersize for plots
    cmdline_args = parser.parse_args()
    if not assume_sources:
        cmdline_args.catalogue = (
            cmdline_args.catalogue_pattern
            %
            get_fname_pattern_substitutions(cmdline_args.frame_fname)
        )
    if cmdline_args.slice is None:
        cmdline_args.slice = [parse_slice('x = 0 +- 0.2'),
                              parse_slice('y = 0 +- 0.2')]
    return cmdline_args
#pylint: enable=too-many-statements


def get_trans_fname(frame_fname, trans_pattern):
    """
    Return the filename of the trans file corresponding to a given frame.

    Args:
        frame_fname(str):    The filename of the FITS frame used for plotting.

        trans_pattern(str):    The pattern for the filename of the
            transformation file specified on the command line.

    Returns:
        str:
            The filename of the transformation file.
    """

    trans_fname = (
        trans_pattern
        %
        get_fname_pattern_substitutions(frame_fname)
    )

    if not os.path.exists(trans_fname):
        raise IOError(f'Transformation file {trans_fname!r} does not exist!')

    return trans_fname


def get_source_positions(catalogue_fname, trans_fname, image_resolution):
    """
    Return [(x, y), ...] positions of the sources in the image.

    Args:
        catalogue_fname(str):    The filename of a catalogue query containing
            all sources in the image (could contain more).

        trans_fname(str):   The filename of the grtrans file for transforming
            tan-projected (xi, eta) to image positions.

        image_resolution((x_res, y_res)):   The resolution of the input image.
            Used to determine which sources project outside the image.

    Returns:
        [(float, float), ...]:
            A list of the (x, y) positions of the sources in the image.
    """

    def in_image(xy_tuple):
        """True iff the given (x, y) tuple is inside the image."""

        #x and y are perfectly valid names
        #pylint: disable=invalid-name
        x, y = xy_tuple
        return (
            0 < x < image_resolution[0] + 1
            and
            0 < y < image_resolution[1] + 1
        )
        #pylint: enable=invalid-name

    with open(trans_fname, encoding='utf-8') as trans:
        for line in trans:
            if line.startswith('# 2MASS:'):
                field_center = tuple(line.split()[2:4])

    get_xi_eta_cmd = subprocess.Popen(
        [
            'grtrans',
            '--input', catalogue_fname,
            '--wcs', f'tan,degrees,ra={field_center[0]},dec={field_center[1]}',
            '--col-radec', '2,3',
            '--col-out', '2,3',
            '--output', '-'
        ],
        stdout=subprocess.PIPE
    )

    get_x_y_cmd = subprocess.Popen(
        [
            'grtrans',
            '--input', '-',
            '--col-xy', '2,3',
            '--input-transformation', trans_fname,
            '--col-out', '2,3',
            '--output', '-'
        ],
        stdin=get_xi_eta_cmd.stdout,
        stdout=subprocess.PIPE
    )
    get_xi_eta_cmd.stdout.close()
    projected_sources = get_x_y_cmd.communicate()[0]

    return list(
        filter(
            in_image,
            (
                tuple(float(v) for v in line.split()[1:3])
                for line in projected_sources.strip().split(b'\n')
            )
        )
    )


def get_source_info(*,
                    pixel_array,
                    stddev_array,
                    mask_array,
                    source_positions,
                    aperture,
                    bg_radii):
    """
    Return field array containing source positions, fluxes and backgrounds.

    Args:
        pixel_array (2-D array like):    The measured values of the image
            pixels.

        stddev_array (2-D array like):    The estimated variance of the image
            pixels.

        mask_array (2-D array like):    Quality flags for the image pixels.

        source_positions:    The return value from get_source_positions().

        aperture:    The size of the aperture to use for measuring the flux.

        bg_radii((float, float)):    The inner and outer radius to use for the
            background annulus.

    Returns:
        (scipy field array):
            All relevant source information in the following fields:

                - x: The x coordinates of the sources.

                - y: The y coordinates of the sources.

                - flux: The measured fluxes of the sources.

                - flux_err: Estimated error of the flux.

                - bg: The measured backgrounds of the sources.

                - bg_err: Estimated error of the background.

                - bg_npix: The number of pixels used in determining the
                  background.
    """

    def add_flux_info(result, measure_background):
        """Measure the flux of the sources and add to result."""

        fit_star_shape = FitStarShape(
            mode='PSF',
            shape_terms='{1}',
            grid=[[-1.1 * aperture, 1.1 * aperture],
                  [-1.1 * aperture, 1.1 * aperture]],
            initial_aperture=aperture,
            bg_min_pix=0,
            src_min_pix=0,
            src_min_signal_to_noise=0.0,
            src_max_aperture=1000.0
        )
        result_tree = fit_star_shape.fit(
            [
                (
                    pixel_array,
                    stddev_array,
                    mask_array,
                    result
                )
            ],
            [measure_background]
        )
        get_flux = SubPixPhot(apertures=[aperture])
        get_flux(
            (pixel_array, stddev_array, mask_array),
            result_tree,
        )
        result['flux'] = flux_from_magnitude(
            result_tree.get(quantity='apphot.mag.0.0',
                            shape=result.shape),
            get_flux.configuration['magnitude_1adu']
        )

    result = scipy.empty(len(source_positions),
                         dtype=[('ID', 'S5'),
                                ('x', scipy.float64),
                                ('y', scipy.float64),
                                ('flux', scipy.float64),
                                ('bg', scipy.float64),
                                ('bg_err', scipy.float64),
                                ('bg_npix', scipy.uint64),
                                ('enabled', scipy.float64)])
    result['enabled'][:] = True
    src_x = scipy.fromiter((pos[0] for pos in source_positions), dtype=c_double)
    src_y = scipy.fromiter((pos[1] for pos in source_positions), dtype=c_double)
    for int_id in range(result.size):
        result['ID'][int_id] = f'{int_id:05d}'
    result['x'] = src_x
    result['y'] = src_y

    measure_background = BackgroundExtractor(pixel_array, *bg_radii)
    result['bg'], result['bg_err'], result['bg_npix'] = measure_background(
        src_x,
        src_y
    )

    add_flux_info(result, measure_background)

    return result


#No clean way to simplify found.
#pylint: disable=too-many-locals
def find_pixel_offsets(sources,
                       prf_range,
                       image_resolution,
                       crowding_distance,
                       plot=False):
    """
    Return the positions of pixels within prf range relative to source centers.

    Args:
        sources:    The return value of get_source_info().

        prf_range:    See `--prf-range` command line argument.

        image_resolution(int, int):    The x and y resolution of the image under
            investigation.

        crowding_distance(float):    The minimum distance between a source and
            its closest neighbor to still consider the source isolated.

        plot(bool):    If True, show plots of the offsets and norm images.

    Returns:
        2-D field array:

            * `x_off`: The field giving the offset of the pixel center from
              the source position in the x direction

            * `y_off`: The field giving the offset of the pixel center from
              the source position in the y direction

            * `norm`: The normalization to use for the pixel response.

            Pixels that are not within range of any PSF or that are within more
            than one PSF's range have `nan` entries.
    """

    result = scipy.full(
        image_resolution,
        scipy.nan,
        dtype=[('x_off', scipy.float64),
               ('y_off', scipy.float64),
               ('norm', scipy.float64),
               ('zero_point', scipy.float64)]
    )
    shared = scipy.full(image_resolution, False, dtype=bool)

    #False positive.
    #pylint: disable=no-member
    source_tree = scipy.spatial.cKDTree(scipy.c_[sources['x'], sources['y']])
    #pylint: enable=no-member
    crowded_flags = source_tree.query_ball_point(source_tree.data,
                                                 crowding_distance,
                                                 return_length=True) > 1

    for this_source, crowded in zip(sources, crowded_flags):
        min_x = int(max(scipy.floor(this_source['x'] - prf_range[2]), 0))
        max_x = int(
            min(
                scipy.ceil(this_source['x'] - prf_range[2] + prf_range[0]),
                image_resolution[1]
            )
        )
        min_y = int(max(scipy.floor(this_source['y'] - prf_range[3]), 0))
        max_y = int(
            min(
                scipy.ceil(this_source['y'] - prf_range[3] + prf_range[1]),
                image_resolution[0]
            )
        )
        result_patch = result[min_y : max_y, min_x : max_x]
        if crowded:
            shared[min_y : max_y, min_x : max_x] = True
        else:
            shared[
                min_y : max_y, min_x : max_x
            ][
                scipy.isfinite(result_patch['x_off'])
            ] = True
        result_patch['x_off'] = (scipy.arange(min_x, max_x)
                                 -
                                 this_source['x']
                                 +
                                 0.5)

        result_patch['y_off'].transpose()[:] = (scipy.arange(min_y, max_y)
                                                -
                                                this_source['y']
                                                +
                                                0.5)
        result_patch['zero_point'] = this_source['bg']
        result_patch['norm'] = this_source['flux']
    result[shared] = scipy.nan
    if plot:
        pyplot.imshow(result['x_off'], origin='lower')
        pyplot.colorbar()
        pyplot.show()

        pyplot.imshow(result['y_off'], origin='lower')
        pyplot.colorbar()
        pyplot.show()

        pyplot.imshow(result['norm'], origin='lower')
        pyplot.colorbar()
        pyplot.show()

    return result
#pylint: enable=too-many-locals


def get_prf_data(pixel_values,
                 pixel_stddev,
                 pixel_offsets,
                 error_threshold=0.1):
    """
    Return the PRF measurements for (a subset of) the image.

    Args:
        pixel_values(2-D float array):    The calibrated pixel responses from
            the image to include in the plot.

        pixel_stddev(2-D float array):    The estimated standard deviation of
            `pixel_values`.

        pixel_offsets:    The slice of the return value of find_pixel_offsets()
            corresponding to `pixel_values`.

        error_threshold(float):    See `--error-threshold` command line
            argument

    Returns:
        (2-D float array, 2-D float array, 2-D float array, 2-D float array):

            * The x-offsets of the points at which PRF measurements are
              available.

            * The y-offsets of the points at which PRF measurements are
              available.

            * The measured normalized PRF at the available offsets

            * estimated errors of the PRF measurements.
    """

    prf_measurements = (
        (
            pixel_values
            -
            pixel_offsets['zero_point']
        )
        /
        pixel_offsets['norm']
    )

    prf_errors = (
        pixel_stddev
        /
        pixel_offsets['norm']
    )
    #False positive
    #pylint: disable=assignment-from-no-return
    include = scipy.logical_and(scipy.isfinite(prf_measurements),
                                scipy.isfinite(prf_errors))
    include = scipy.logical_and(include, prf_errors < error_threshold)
    #pylint: enable=assignment-from-no-return
    return scipy.stack((
        pixel_offsets['x_off'][include],
        pixel_offsets['y_off'][include],
        prf_measurements[include],
        prf_errors[include]
    ))


#TODO: see if can be fixed
#pylint: disable=too-many-locals
def plot_prf_slice(prf_data,
                   spline,
                   *,
                   label,
                   x_offset=None,
                   y_offset=None,
                   thickness=0.1,
                   error_scale=0.1,
                   points_color='k',
                   marker_size=2,
                   **binning):
    """
    Plot a slice of the PRF.

    Args:
        prf_data:    The return value of get_prf_data().

        x_offset/y_offset(float):    Plot a slice of constant x or y (determined
            by the argument name) offset from the source center.

        thickness(float):    Points with x or y within `thickness` of the
            specified offset are included in the plot.

        error_scale(float):    See `--error-scale` command line argument.

        error_threshold(float):    See `--error-threshold` command line
            argument.

    Returns:
        None
    """

    assert (
        (x_offset is None and y_offset is not None)
        or
        (x_offset is not None and y_offset is None)
    )

    assert len(prf_data) == 4
    for entry in prf_data[1:]:
        assert entry.shape == prf_data[0].shape

    if x_offset is None:
        plot_pixel_indices = scipy.nonzero(
            scipy.fabs(prf_data[1] - y_offset) < thickness
        )
        if spline is not None:
            spline_x = scipy.linspace(prf_data[0].min(), prf_data[0].max(), 300)
            spline_y = spline(spline_x, y_offset).flatten()
    else:
        plot_pixel_indices = scipy.nonzero(
            scipy.fabs(prf_data[0] - x_offset) < thickness
        )
        if spline is not None:
            spline_x = scipy.linspace(prf_data[1].min(), prf_data[1].max(), 300)
            spline_y = spline(x_offset, spline_x).flatten()

    plot_x = prf_data[
        0 if x_offset is None else 1
    ][
        plot_pixel_indices
    ]

    plot_y = prf_data[2][plot_pixel_indices]

    plot_err_y = prf_data[3][plot_pixel_indices]

    if plot_x.size == 0:
        return

    pyplot.errorbar(plot_x,
                    plot_y,
                    plot_err_y * error_scale,
                    fmt='o',
                    color=points_color,
                    markersize=1.5 * marker_size,
                    elinewidth=(marker_size / 2),
                    zorder=10,
                    label=label)
    if spline is not None:
        pyplot.plot(spline_x,
                    spline_y,
                    '-',
                    color='black',
                    linewidth=marker_size,
                    zorder=20,
                    alpha=0.85)

    if binning:
        pyplot.plot(
            scipy.stats.binned_statistic(plot_x,
                                         plot_x,
                                         **binning)[0],
            scipy.stats.binned_statistic(plot_x,
                                         plot_y,
                                         **binning)[0],
            'o',
            markerfacecolor=points_color,
            markeredgecolor='black',
            markersize=(5 * marker_size),
            linewidth=3,
            zorder=30,
            alpha=0.7
        )
#pylint: enable=too-many-locals


def plot_entire_prf(cmdline_args,
                    image_slices,
                    grid_x,
                    grid_y,
                    sources=None):
    """
    Plots the entire PRF on a 3d axes and displays it to the user.

    Args:
        cmdline_args:    The parsed command line arguments.

        image_slices:    How to split the image when plotting (the return
            value of get_image_slices())

        grid_x:     The PSF x-grid, pulled from star_shape_grid

        grid_y:     The PSF y-grid, pulled from star_shape_grid

        sources(scipy.array or None):    If not None, it should be an array with
            equivalent structure to get_source_info(). If None,
            get_source_info() is used to initialize it.


    Returns:
        None
    """

    def output_plot(x_image_slice, y_image_slice, x_index, y_index):
        """Save or display the currently set up plot."""

        if cmdline_args.save_prf_plot is None:
            pyplot.show()
            return

        substitutions = get_fname_pattern_substitutions(
            cmdline_args.frame_fname
        )
        substitutions['x_split'] = str(x_image_slice) + '_' + str(x_index)
        substitutions['y_split'] = str(y_image_slice) + '_' + str(y_index)
        fname = cmdline_args.save_prf_plot % substitutions
        file_exists = prepare_file_output(fname,
                                          cmdline_args.overwrite_plots,
                                          True,
                                          True)

        if cmdline_args.skip_existing_plots and file_exists:
            return

        pyplot.savefig(fname,
                       dpi=cmdline_args.figure_dpi)
        pyplot.cla()

    def plot_slice(frame, pixel_offsets, x_image_slice, y_image_slice):
        """Plot the PRF of one slice of the image."""

        print('first_hdu='+repr(first_hdu))
        print('x_image_slice='+repr(x_image_slice))
        print('y_image_slice='+repr(y_image_slice))
        print('frame_length='+repr(len(frame)))

        prf_data = get_prf_data(
            frame[first_hdu].data[y_image_slice, x_image_slice],
            frame[first_hdu + 1].data[y_image_slice, x_image_slice],
            pixel_offsets[y_image_slice, x_image_slice],
            cmdline_args.error_threshold
        )
        plot_x = prf_data[0]
        plot_y = prf_data[1]
        plot_z = prf_data[2]

        # 3D
        # minz = numpy.amin(plot_z)
        # maxz = numpy.amax(plot_z)
        # norm_z = (plot_z - plot_z.min()) / plot_z.max()
        # colors = cm.hsv(norm_z)
        # ax = pyplot.axes(projection='3d')
        # ax.scatter3D(plot_x,
        #              plot_y,
        #              plot_z,
        #              cmap=colors,
        #              c=colors,
        #              s=10,
        #              marker='o',
        #              alpha=0.3,
        #              vmin=minz,
        #              vmax=maxz
        #              )

        pyplot.scatter(plot_x,
                       plot_y,
                       cmap='hsv',
                       c=plot_z,
                       s=10,
                       marker='o',
                       alpha=0.3,
                       vmin=-0.001,
                       vmax=0.15
                       )
        pyplot.colorbar()
        pyplot.xticks(grid_x)
        pyplot.yticks(grid_y)
        pyplot.grid(grid_y, linestyle='--', color='k')
        pyplot.grid(grid_x, linestyle='--', color='k')

    with fits.open(cmdline_args.frame_fname, 'readonly') as frame:
        # False positive
        # pylint: disable=no-member
        if frame[0].header['NAXIS']:
            image_resolution = (frame[0].header['NAXIS2'],
                                frame[0].header['NAXIS1'])
            first_hdu = 0
        else:
            image_resolution = (frame[1].header['NAXIS2'],
                                frame[1].header['NAXIS1'])
            first_hdu = 1
        # pylint: enable=no-member

        assert image_resolution == frame[first_hdu].data.shape

        if sources is None:
            sources = get_source_info(
                pixel_array=frame[first_hdu].data.astype(float),
                stddev_array=frame[first_hdu + 1].data.astype(float),
                mask_array=frame[first_hdu + 2].data.astype(c_char),
                source_positions=get_source_positions(
                    cmdline_args.catalogue,
                    get_trans_fname(cmdline_args.frame_fname,
                                    cmdline_args.trans_pattern),
                    image_resolution
                ),
                aperture=cmdline_args.flux_aperture,
                bg_radii=cmdline_args.background_annulus
            )

        pixel_offsets = find_pixel_offsets(sources,
                                           cmdline_args.prf_range,
                                           image_resolution,
                                           2.0 * cmdline_args.flux_aperture)

        for x_image_slice, y_image_slice, x_index, y_index in image_slices:

            plot_slice(frame, pixel_offsets, x_image_slice, y_image_slice)
            output_plot(x_image_slice,
                        y_image_slice,
                        x_index,
                        y_index)

            #TODO need to manually set vmin and vmax make color have a
            #wider range try to find color limit, should make them quantiles
            #functions

def plot_3d_prf(cmdline_args, meshgrid_x, meshgrid_y, prf_slice_splines):
    """
    Plot a 3D slice of the PRF.

    Args:
        cmdline_args: The parse command line configuration.

        meshgrid_x, meshgrid_y:   The meshgrid generated for evaluating the PRF
            made from the grid points of the star_shape_grid

        prf_slice_splines:    The evaluated PRF at the meshgrid points extracted
            from slice_splines

    Returns:
        None
    """

    def plot_3d_slice(spline, meshgrid_x, meshgrid_y, stride, linewidths):
        """Plot the 3D PRF of one slice of the image."""
        pyplot.axes(
            projection='3d'
        ).contour3D(
            meshgrid_x,
            meshgrid_y,
            spline.reshape(meshgrid_x.shape),
            stride,
            linewidths=linewidths
        )

    def get_3d_plot_tasks(cmdline_args):
        """
        Generator for all the 3D plots that must be created.

        Args:
            cmdline_args:    The parse command line configuration.

        Yields:
            [(str on None), ...]: List of the filenames to save the plot if each
                slice is saved separately. Appropriate length of Nones
                otherwise.

            str or None: The filename to save a combined plot of all pieces of
                the image together or None if individual images are saved
                separately.
        """

        def get_3d_image_slice_fname(x_image_slice,
                                     y_image_slice,
                                     x_index,
                                     y_index,
                                     fname_substitutions):
            """Return a single ently of the plotting tasks per image silce."""

            fname_substitutions['XSPLIT'] = (str(x_image_slice) + str(x_index)
                                             if x_image_slice else
                                             None)
            fname_substitutions['YSPLIT'] = (str(y_image_slice) + str(y_index)
                                             if y_image_slice else
                                             None)
            return (
                cmdline_args.save_3d_spline_plot % fname_substitutions
                if cmdline_args.plot_multi_image else
                None
            )


        fname_substitutions = get_fname_pattern_substitutions(
            cmdline_args.frame_fname
        )
        image_slice_list = get_image_slices(cmdline_args.split_image,
                                            cmdline_args.discard_image_boundary)
        yield (
            [get_3d_image_slice_fname(*image_slice, fname_substitutions)
             for image_slice in image_slice_list],
            (
                None if cmdline_args.plot_multi_image
                else cmdline_args.save_3d_spline_plot % fname_substitutions
            )
        )

    allow_existing_plot = (cmdline_args.overwrite_plots
                           or
                           cmdline_args.skip_existing_plots)

    for slice_fnames, combined_fname in get_3d_plot_tasks(
            cmdline_args
    ):
        if (
                combined_fname is not None
                and
                prepare_file_output(combined_fname,
                                    allow_existing_plot,
                                    True,
                                    cmdline_args.overwrite_plots)
        ):
            continue

        if cmdline_args.plot_y_range is not None:
            pyplot.ylim(*cmdline_args.plot_y_range)

        for (
                spline,
                individual_fname
        ) in zip(
            prf_slice_splines,
            slice_fnames
        ):
            if (
                    individual_fname is not None
                    and
                    prepare_file_output(individual_fname,
                                        allow_existing_plot,
                                        True,
                                        cmdline_args.overwrite_plots)
            ):
                continue

            plot_3d_slice(
                spline,
                meshgrid_x,
                meshgrid_y,
                cmdline_args.plot_3d_spline_contour_stride,
                linewidths=cmdline_args.plot_3d_spline_contour_linewidths
            )

            if cmdline_args.plot_multi_image:

                if individual_fname is None:
                    pyplot.show()
                else:
                    pyplot.savefig(individual_fname,
                                   dpi=cmdline_args.figure_dpi)
                    pyplot.cla()

        if cmdline_args.plot_multi_image:
            continue

        pyplot.legend()

        if combined_fname is None:
            pyplot.show()
        else:
            pyplot.savefig(combined_fname,
                           dpi=cmdline_args.figure_dpi)
            pyplot.cla()


def get_image_slices(splits, inner_only):
    """
    Return [(x-slice, y-slice, x_index, y_index), ...] per `--split-image` arg.

    Args:
        splits:    The parsed image splits directly from the command line.

        inner_only(bool):    If True, slices in contact with the outer border
            of the image are discarded.

    Returns:
        [(x-slice, y-slice, x_index, y_index), ...]:
            Array slices directly useable to index the image implementing the
            splitting specified on the command line.
    """

    split_slices = {'x': [slice(0, None)], 'y': [slice(0, None)]}

    for direction, value in sorted(splits):
        split_slices[direction][-1] = slice(split_slices[direction][-1].start,
                                            value)
        split_slices[direction].append(slice(value, None))

    if inner_only:
        for direction in 'xy':
            split_slices[direction] = split_slices[direction][1:-1]

    print('Split slices: ' + repr(split_slices))

    return [
        (x_split, y_split, x_index, y_index)
        for x_index, x_split in enumerate(split_slices['x'])
        for y_index, y_split in enumerate(split_slices['y'])
    ]


class AlglibSpline:
    """Wrap alglib based splines so they can be evaluated on scipy arrays."""

    def __init__(self, prf_data, resolution, penalty, domain):
        """
        Fit the spline on the given data.

        Args:
            prf_data:    The return value of get_prf_data().

            resolution(int, int):    The grid resolution to use when building
                the spline.

            penalty(float):    The nonlinearity penalty when fitting the spline.

            domain(float, float, float, float):    The area (xmin, xmax, ymin,
                ymax) over which the spline will be derived.

        Returns:
            None
        """

        builder = xalglib.spline2dbuildercreate(1)
        xalglib.spline2dbuildersetpoints(builder,
                                         prf_data[:3].T.tolist(),
                                         prf_data[0].size)
        xalglib.spline2dbuildersetgrid(builder, *resolution)
        xalglib.spline2dbuildersetarea(builder, *domain)
        xalglib.spline2dbuildersetalgoblocklls(builder, penalty)

        spline = xalglib.spline2dfit(builder)[0]

        ny_nodes, nx_nodes, nvals, coef_table = xalglib.spline2dunpackv(spline)

        assert nvals == 1

        self.x_nodes = (
            [coef_table[i][0] for i in range(nx_nodes - 1)]
            +
            [coef_table[nx_nodes - 2][1]]
        )
        self.y_nodes = (
            [coef_table[j * (nx_nodes - 1)][2] for j in range(ny_nodes - 1)]
            +
            [coef_table[(ny_nodes - 2) * (nx_nodes - 1)][3]]
        )

        self.spline_eval = scipy.vectorize(
            lambda x, y: xalglib.spline2dcalc(spline, x, y)
        )

    def __call__(self, x, y):
        """Evaluate the spline."""

        return self.spline_eval(x, y)

    def plot_grid_boundaries(self, direction, color):
        """
        Add vertical lines to the current axis showing the grid along direction.

        Args:
            direction(str):    Either `'x'`, or `'y'`, selecting which direction
                to show the grid lines along.

        Returns:
            None
        """

        for node in getattr(self, direction + '_nodes'):
            pyplot.axvline(x=node, color=color)


#TODO: simplify later
#pylint: disable=too-many-locals
def pad_prf_data(prf_data, cmdline_args):
    """Return PRF data zero-padded for fitting and the spline fit domain."""

    #False positive
    #pylint: disable=assignment-from-no-return
    keep_prf_data = scipy.logical_and(
        scipy.logical_and(
            prf_data[0] > -cmdline_args.prf_range[2],
            prf_data[0] < cmdline_args.prf_range[0] - cmdline_args.prf_range[2]
        ),
        scipy.logical_and(
            prf_data[1] > -cmdline_args.prf_range[3],
            prf_data[1] < cmdline_args.prf_range[1] - cmdline_args.prf_range[3]
        )
    )
    #pylint: enable=assignment-from-no-return
    num_keep = keep_prf_data.sum()

    x_padding = cmdline_args.prf_range[0] * cmdline_args.spline_pad_fraction
    y_padding = cmdline_args.prf_range[1] * cmdline_args.spline_pad_fraction
    domain = (
        -cmdline_args.prf_range[2] - x_padding,
        cmdline_args.prf_range[0] - cmdline_args.prf_range[2] + x_padding,
        -cmdline_args.prf_range[3] - y_padding,
        cmdline_args.prf_range[1]-cmdline_args.prf_range[3] + y_padding
    )
    middle_npoints = int(cmdline_args.spline_pad_npoints
                         /
                         cmdline_args.spline_pad_fraction)
    padded_prf_data = scipy.empty(
        (
            4,
            (
                num_keep
                +
                4
                *
                cmdline_args.spline_pad_npoints
                *
                (
                    middle_npoints
                    +
                    cmdline_args.spline_pad_npoints
                )
            )
        ),
        dtype=float
    )
    padded_prf_data[:, : num_keep] = prf_data[:, keep_prf_data]
    padded_prf_data[2, num_keep : ] = 0.0
    padded_prf_data[3, num_keep : ] = (
        prf_data[3, keep_prf_data].min()
        /
        num_keep
    )
    xy_padding_start = num_keep
    corner_x, corner_y = scipy.meshgrid(
        scipy.linspace(0, x_padding, cmdline_args.spline_pad_npoints),
        scipy.linspace(0, y_padding, cmdline_args.spline_pad_npoints)
    )

    mid_x_x, mid_x_y = scipy.meshgrid(
        scipy.linspace(
            -cmdline_args.prf_range[2],
            cmdline_args.prf_range[0] - cmdline_args.prf_range[2],
            middle_npoints + 2
        )[1:-1],
        scipy.linspace(
            0,
            y_padding,
            cmdline_args.spline_pad_npoints
        )
    )

    mid_y_x, mid_y_y = scipy.meshgrid(
        scipy.linspace(
            0,
            x_padding,
            cmdline_args.spline_pad_npoints
        ),
        scipy.linspace(
            -cmdline_args.prf_range[3],
            cmdline_args.prf_range[1] - cmdline_args.prf_range[3],
            middle_npoints + 2
        )[1:-1]
    )

    corner_x = corner_x.flatten()
    corner_y = corner_y.flatten()
    mid_x_x = mid_x_x.flatten()
    mid_x_y = mid_x_y.flatten()
    mid_y_x = mid_y_x.flatten()
    mid_y_y = mid_y_y.flatten()

    add_mid_x = True
    for x_offset in [domain[0], domain[1] - x_padding]:
        padded_prf_data[
            0,
            xy_padding_start: xy_padding_start + mid_y_x.size
        ] = mid_y_x + x_offset
        padded_prf_data[
            1,
            xy_padding_start: xy_padding_start + mid_y_x.size
        ] = mid_y_y
        xy_padding_start += mid_y_x.size
        for y_offset in [domain[2], domain[3] - y_padding]:
            if add_mid_x:
                padded_prf_data[
                    0,
                    xy_padding_start: xy_padding_start + mid_x_x.size
                ] = mid_x_x
                padded_prf_data[
                    1,
                    xy_padding_start: xy_padding_start + mid_x_x.size
                ] = mid_x_y + y_offset
                xy_padding_start += mid_x_x.size
            padded_prf_data[
                0,
                xy_padding_start: xy_padding_start + corner_x.size
            ] = corner_x + x_offset
            padded_prf_data[
                1,
                xy_padding_start: xy_padding_start + corner_x.size
            ] = corner_y + y_offset
            xy_padding_start += corner_x.size
        add_mid_x = False

    return padded_prf_data, domain
#pylint: enable=too-many-locals


def fit_spline(prf_data, domain, cmdline_args):
    """Return the best-fit spline to the PRF per the command line config."""

    if cmdline_args.spline_method == 'none':
        return None

    if cmdline_args.spline_method == 'scipy':
        return SmoothBivariateSpline(
            prf_data[0],
            prf_data[1],
            prf_data[2],
            1.0 / prf_data[3],
            s=cmdline_args.spline_smoothing * prf_data[3].size,
            bbox=domain
        )

    return AlglibSpline(
        prf_data,
        cmdline_args.spline_resolution,
        cmdline_args.spline_smoothing,
        domain=domain
    )

def get_plot_tasks(cmdline_args):
    """
    Generator for all the plots that must be created.

    Args:
        cmdline_args:    The parse command line configuration.

    Yields:
        dict: The parsed PRF slice directly taken from the command line to plot

        [(str on None), ...]: List of the filenames to save the plot if each
            slice is saved separately. Appropriate length of Nones otherwise.

        str or None: The filanem to save a combined plot of all pieces of the
           image together or None if individual images are saved separately.
    """

    def get_image_slice_fname(x_image_slice,
                              y_image_slice,
                              x_index,
                              y_index,
                              fname_substitutions):
        """Return a single ently of the plotting tasks per image silce."""

        fname_substitutions['XSPLIT'] = (str(x_image_slice) + str(x_index)
                                         if x_image_slice else
                                         None)
        fname_substitutions['YSPLIT'] = (str(y_image_slice) + str(y_index)
                                         if y_image_slice else
                                         None)
        return (
            cmdline_args.save_plot_pattern % fname_substitutions
            if cmdline_args.plot_multi_image else
            None
        )

    fname_substitutions = get_fname_pattern_substitutions(
        cmdline_args.frame_fname
    )
    image_slice_list = get_image_slices(cmdline_args.split_image,
                                        cmdline_args.discard_image_boundary)

    for plot_slice in cmdline_args.slice:
        slice_direction = ('x' if 'x_offset' in plot_slice else 'y')
        fname_substitutions['SLICEDIR'] = slice_direction
        fname_substitutions['SLICEOFFSET'] = plot_slice[slice_direction
                                                        +
                                                        '_offset']
        yield (
            plot_slice,
            [get_image_slice_fname(*image_slice, fname_substitutions)
             for image_slice in image_slice_list],
            (
                None if cmdline_args.plot_multi_image
                else cmdline_args.save_plot_pattern % fname_substitutions
            )
        )

def list_plot_filenames(cmdline_args):
    """List the filenames of all the plots that will be generated."""

    assert cmdline_args.save_plot_pattern is not None
    if not cmdline_args.plot_multi_image:
        return [task[1] for task in get_plot_tasks(cmdline_args)]

    result = []
    for _, slice_fnames, outer_plot_fname in get_plot_tasks(cmdline_args):
        assert outer_plot_fname is None
        result.extend(slice_fnames)

    return result

#TODO: see if it can be simplified
#pylint: disable=too-many-branches
def show_plots(slice_prf_data, slice_splines, cmdline_args, append=False):
    """
    Generate the plots and display them to the user.

    Args:
        slice_prf_data(iterable):    List of the result of either get_prf_data()
            or pad_prf_data() for each image slice.

        slice_splines(iterable):    List of best-fit splines to the data for
            each image slice. Must support evaluation using arrays of x & y
            positions.

        cmdline_args:    The parsed command line arguments.

    Returns
        None
    """

    allow_existing_plot = (cmdline_args.overwrite_plots
                           or
                           cmdline_args.skip_existing_plots)
    for plot_slice, slice_fnames, combined_fname in get_plot_tasks(
            cmdline_args
    ):
        print('plot_slice: ' + repr(plot_slice))
        print('Slice filenames: ' + repr(slice_fnames))
        print('combined_fname: ' + repr(combined_fname))
        if (
                combined_fname is not None
                and
                prepare_file_output(combined_fname,
                                    allow_existing_plot,
                                    True,
                                    cmdline_args.overwrite_plots)
        ):
            continue

        if cmdline_args.plot_y_range is not None:
            pyplot.ylim(*cmdline_args.plot_y_range)

        if len(slice_fnames) == 1 and slice_fnames[0] is None:
            slice_fnames *= len(slice_prf_data)
        for (
                (prf_data, label),
                spline,
                individual_fname,
                color
        ) in zip(
            slice_prf_data,
            slice_splines,
            slice_fnames,
            kelly_colors[2:]
        ):
            if (
                    individual_fname is not None
                    and
                    prepare_file_output(individual_fname,
                                        allow_existing_plot,
                                        True,
                                        cmdline_args.overwrite_plots)
            ):
                continue

            plot_prf_slice(
                prf_data[0],
                spline,
                label=(label if not append else None),
                error_scale=cmdline_args.error_scale,
                points_color=color,
                marker_size=cmdline_args.marker_size,
                **plot_slice,
                **cmdline_args.add_binned
            )
            try:
                spline.plot_grid_boundaries(
                    ('x' if 'y_offset' in plot_slice else 'y'),
                    color=color
                )
            except AttributeError:
                pass
            if cmdline_args.plot_multi_image:
                assert not append
                pyplot.axhline(y=0)
                pyplot.legend()
                pyplot.xlabel('pixel center - source center [pix]')
                pyplot.ylabel('normalized pixel response')

                if individual_fname is None:
                    pyplot.show()
                else:
                    pyplot.savefig(individual_fname,
                                   dpi=cmdline_args.figure_dpi)
                    pyplot.cla()

        if cmdline_args.plot_multi_image or append:
            continue

        print('Saving: ' + repr(combined_fname))
        pyplot.axhline(y=0)
        pyplot.xlabel('pixel center - source center [pix]')
        pyplot.ylabel('normalized pixel response')
        pyplot.legend()

        if combined_fname is None:
            pyplot.show()
        else:
            pyplot.savefig(combined_fname,
                           dpi=cmdline_args.figure_dpi)
            pyplot.cla()
#pylint: enable=too-many-branches


def extract_pixel_data(cmdline_args, image_slices, sources=None):
    """
    Get the pixel level data from the input image required for the plot.

    Args:
        cmdline_args:    The parsed command line arguments.

        image_slices:    How to split the image when plotting (the return value
            of get_image_slices()).

        sources(scipy.array or None):    If not None, it should be an array with
            equivalent structure to get_source_info(). If None,
            get_source_info() is used to initialize it.

    Returns:
        List:
            A list of the PRF data (see return of get_prf_data()) for each
            image slice each entry contains the data and a plot label.
    """

    with fits.open(cmdline_args.frame_fname, 'readonly') as frame:
        #False positive
        #pylint: disable=no-member
        if frame[0].header['NAXIS']:
            image_resolution = (frame[0].header['NAXIS2'],
                                frame[0].header['NAXIS1'])
            first_hdu = 0
        else:
            image_resolution = (frame[1].header['NAXIS2'],
                                frame[1].header['NAXIS1'])
            first_hdu = 1

        assert image_resolution == frame[first_hdu].data.shape

        if sources is None:
            trans_fname = get_trans_fname(cmdline_args.frame_fname,
                                          cmdline_args.trans_pattern)
            #pylint: enable=no-member
            sources = get_source_info(
                pixel_array=frame[first_hdu].data.astype(float),
                stddev_array=frame[first_hdu + 1].data.astype(float),
                mask_array=frame[first_hdu + 2].data.astype(c_char),
                source_positions=get_source_positions(cmdline_args.catalogue,
                                                      trans_fname,
                                                      image_resolution),
                aperture=cmdline_args.flux_aperture,
                bg_radii=cmdline_args.background_annulus
            )
            #pylint: disable=no-member

        pixel_offsets = find_pixel_offsets(sources,
                                           cmdline_args.prf_range,
                                           image_resolution,
                                           2.0 * cmdline_args.flux_aperture)
        return [
            (
                pad_prf_data(
                    get_prf_data(
                        frame[first_hdu].data[y_image_slice, x_image_slice],
                        frame[first_hdu + 1].data[y_image_slice, x_image_slice],
                        pixel_offsets[y_image_slice, x_image_slice],
                        cmdline_args.error_threshold,
                    ),
                    cmdline_args
                ),
                (
                    #This is more readable
                    #pylint: disable=consider-using-f-string
                    '(%f, %f)'
                    %
                    (
                        (
                            x_image_slice.start
                            +
                            (x_image_slice.stop or image_resolution[1])
                        ) / 2.0,
                        (
                            y_image_slice.start
                            +
                            (y_image_slice.stop or image_resolution[0])
                        ) / 2.0
                    )
                    #pylint: enable=consider-using-f-string
                )
            )
            for x_image_slice, y_image_slice, x_index, y_index in image_slices
        ]


def main(cmdline_args):
    """Avoid polluting global namespace."""

    image_slices = get_image_slices(cmdline_args.split_image,
                                    cmdline_args.discard_image_boundary)
    slice_prf_data = extract_pixel_data(cmdline_args, image_slices)
    slice_splines = [
        fit_spline(prf_data, domain, cmdline_args)
        for prf_data, domain in slice_prf_data
    ]

    show_plots(slice_prf_data,
               slice_splines,
               cmdline_args)

    #TODO fix this to, broken atm fix sources

    # if cmdline_args.plot_entire_prf:
    #     plot_entire_prf(cmdline_args,
    #                     image_slices,
    #                     sources=None)


if __name__ == '__main__':
    main(parse_command_line())
