#!/usr/bin/env python3

"""Create some plots helpful for picking PSF grid."""

from argparse import ArgumentParser
import os.path
import subprocess
from ctypes import c_double, c_char

import scipy
import scipy.spatial
import scipy.stats
from matplotlib import pyplot
from astropy.io import fits

from superphot import BackgroundExtractor, FitStarShape, SubPixPhot

def parse_command_line():
    """Return the command line arguments as attributes of an object."""

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

    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        'frame_fname',
        help='The full path of the FITS file to use for creating the plots.'
    )
    parser.add_argument(
        '--trans-pattern', '-t',
        default=os.path.join('%(frame_dir)s',
                             '..',
                             'ASTROM',
                             '%(frame_id)s.trans'),
        help="A pattern with substitutions `'%(frame_dir)s'` (directory "
        "containing the frame), and `'%(frame_id)s'` (base filename of the "
        "frame without the `fits` or `fits.fz` extension) that expands to "
        "the `.trans` file corresponding to the input frame. Default: "
        "'%(defalut)s'."
    )
    parser.add_argument(
        '--catalogue', '-c',
        default='/data/ETS_HATP32/ASTROM/catalogue.ucac4',
        help='The full path of the catalogue containing all sources in the '
        'image. Default: \'%(default)s\'.'
    )
    parser.add_argument(
        '--prf-range', '-r',
        default=(11.0, 8.0, 4.0, 4.0),
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
        default=[parse_slice('x = 0 +- 0.1'),
                 parse_slice('y = 0 +- 0.1')],
        help='Add more slices to show. Each slice is specified as an offset '
        'along one if the demensions (x or y) and a range around that to '
        'include in the plot. White space around tokens is allowed and '
        'ignored. Example: "x = 0 +- 0.1". Default: %(default)s.'
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
        '--error-scale',
        type=float,
        default=0.1,
        help='Scale the error bars by this value when plotting to make a more '
        'readable plot. Default: %(defalut)s.'
    )
    parser.add_argument(
        '--error_threshold',
        type=float,
        default=0.1,
        help='Points with error bars larger than this are not included in the '
        'plot. Intended to avoid introducing points that are too noisy. '
        'Default: %(default)s'
    )

    return parser.parse_args()

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

    if frame_fname.lower().endswith('.fz'):
        frame_fname = frame_fname[:-len('.fz')]

    assert frame_fname.lower().endswith('.fits')
    frame_fname = frame_fname[:-len('.fits')]
    frame_dir, frame_id = os.path.split(frame_fname)
    trans_fname = trans_pattern % dict(frame_dir=frame_dir,
                                       frame_id=frame_id)
    assert os.path.exists(trans_fname)
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

    with open(trans_fname) as trans:
        for line in trans:
            if line.startswith('# 2MASS:'):
                field_center = tuple(line.split()[2:4])

    get_xi_eta_cmd = subprocess.Popen(
        [
            'grtrans',
            '--input', catalogue_fname,
            '--wcs', 'tan,degrees,ra=%s,dec=%s' % field_center,
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
        print('Created zero PSF fitter')

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
        print('Finished zero PSF mock fit.')
        print('Mask array: ' + repr(mask_array))

        get_flux = SubPixPhot(apertures=[aperture])
        get_flux(
            (pixel_array, stddev_array, mask_array),
            result_tree,
        )
        result['flux'] = 10.0**(
            (
                get_flux.configuration['magnitude_1adu']
                -
                result_tree.get(quantity='apphot.mag.0.0',
                                shape=result.shape)
            ) / 2.5
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
        result['ID'][int_id] = '%5.5d' % int_id
    result['x'] = src_x
    result['y'] = src_y

    measure_background = BackgroundExtractor(pixel_array, *bg_radii)
    result['bg'], result['bg_err'], result['bg_npix'] = measure_background(
        src_x,
        src_y
    )

    print('Adding flux to sources: ' + repr(result))

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

def plot_prf_slice(pixel_values,
                   pixel_stddev,
                   pixel_offsets,
                   *,
                   x_offset=None,
                   y_offset=None,
                   thickness=0.1,
                   error_scale=0.1,
                   error_threshold=0.1,
                   points_color='k'):
    """
    Plot a slice of the PRF.

    Args:
        pixel_values(2-D float array):    The calibrated pixel responses from
            the image to include in the plot.

        pixel_stddev(2-D float array):    The estimated standard deviation of
            `pixel_values`.

        pixel_offsets:    The slice of the return value of find_pixel_offsets()
            corresponding to `pixel_values`.

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

    assert pixel_values.shape == pixel_stddev.shape
    assert pixel_values.shape == pixel_offsets.shape

    if x_offset is None:
        plot_pixel_indices = scipy.nonzero(
            scipy.fabs(pixel_offsets['y_off'] - y_offset) < thickness
        )
    else:
        plot_pixel_indices = scipy.nonzero(
            scipy.fabs(pixel_offsets['x_off'] - x_offset) < thickness
        )

    plot_x = pixel_offsets[
        'x_off' if x_offset is None else 'y_off'
    ][
        plot_pixel_indices
    ]

    plot_y = (
        (
            pixel_values[plot_pixel_indices]
            -
            pixel_offsets['zero_point'][plot_pixel_indices]
        )
        /
        pixel_offsets['norm'][plot_pixel_indices]
    )

    plot_err_y = (pixel_stddev[plot_pixel_indices]
                  /
                  pixel_offsets['norm'][plot_pixel_indices])

    #False positive
    #pylint: disable=assignment-from-no-return
    finite = scipy.logical_and(scipy.isfinite(plot_x),
                               scipy.isfinite(plot_y))
    finite = scipy.logical_and(finite, plot_err_y < error_threshold)
    #pylint: enable=assignment-from-no-return

    plot_x = plot_x[finite]
    plot_y = plot_y[finite]
    plot_err_y = plot_err_y[finite]

    binned_x = scipy.stats.binned_statistic(plot_x,
                                            plot_x,
                                            statistic='median',
                                            bins=50)[0]
    binned_y = scipy.stats.binned_statistic(plot_x,
                                            plot_y,
                                            statistic='median',
                                            bins=50)[0]

    pyplot.errorbar(plot_x,
                    plot_y,
                    plot_err_y * error_scale,
                    fmt='o' + points_color)
    pyplot.plot(binned_x, binned_y, '-ok', markersize=10, linewidth=3)

def get_image_slices(splits):
    """Return a list of (x-slice, y-slice) tuples per `--split-image` arg."""

    split_slices = dict(x=[slice(0, None)], y=[slice(0, None)])

    for direction, value in sorted(splits):
        split_slices[direction][-1] = slice(split_slices[direction][-1].start,
                                            value)
        split_slices[direction].append(slice(value, None))
    return [
        (x_split, y_split)
        for x_split in split_slices['x']
        for y_split in split_slices['y']
    ]

def main():
    """Avoid polluting global namespace."""

    cmdline_args = parse_command_line()
    trans_fname = get_trans_fname(cmdline_args.frame_fname,
                                  cmdline_args.trans_pattern)
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

        #pylint: enable=no-member
        sources = get_source_info(
            pixel_array=frame[first_hdu].data,
            stddev_array=frame[first_hdu + 1].data,
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


        for plot_slice in cmdline_args.slice:
            for (x_image_slice, y_image_slice), color in zip(
                    get_image_slices(
                        cmdline_args.split_image
                    ),
                    'rgb'
            ):
                plot_prf_slice(
                    frame[first_hdu].data[y_image_slice, x_image_slice],
                    frame[first_hdu + 1].data[y_image_slice, x_image_slice],
                    pixel_offsets[y_image_slice, x_image_slice],
                    error_scale=cmdline_args.error_scale,
                    error_threshold=cmdline_args.error_threshold,
                    points_color=color,
                    **plot_slice
                )
            pyplot.axhline(y=0)
            pyplot.show()

if __name__ == '__main__':
    main()
