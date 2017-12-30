"""A collection of functions used by the fitpsf unit tests."""

import sys
sys.path.append('../../PythonPackage')

from fake_image.Image import Image

import pyfits

from math import ceil

def make_image_and_source_list(sources,
                               image_fname,
                               source_list_fname,
                               psf_fit_fname,
                               extra_variables,
                               subpix_map) :
    """
    Create an image and a list of the sources in it ready for psf fitting.

    Args:
        - sources:
            A list of dictionaries with at least the following keywords:
                - x: 
                    The x coordinate of the source center.
                - y:
                    The y coordinate of the source center.
                - psf:
                    An instance of some sub-class of PSFBase giving the
                    sources's PSF. It should already be scaled to the desired
                    flux.
            Additional keywords may be added to the source list and hence
            available as variables for PSF fitting by listing the names in the
            extra_variables argument.

        - image_fname:
            The filename to save the image under. If a file with this name
            exists it is overwritten.

        - source_list_fname:
            The filename to add the source list to. If a file with this name
            exists it appended to.

        - psf_fit_fname:
            The filename to use for PSF fitting results.

        - extra_variables:
            A list of additional keywords from sources to add to the source list
            and the order in which those should be added. The corresponding
            entries in sources must be floating point values.

        - subpix_map:
            The sub-pixel map to impose on the image. For more details see
            same name argument of Image.add_source.

    Returns: None
    """

    min_x = min(s['x'] for s in sources)
    max_x = max(s['x'] for s in sources)
    min_y = min(s['y'] for s in sources)
    max_y = max(s['y'] for s in sources)

    image = Image(int(ceil(max_x + min_x)),
                  int(ceil(max_y + min_y)),
                  background = 1.0)

    src_list = open(source_list_fname, 'a')
    src_list.write('[' + image_fname + ' ' + psf_fit_fname + ']\n')

    src_list_vars = ['x', 'y'] + extra_variables

    for src_id, src in enumerate(sources) :
        image.add_source(x = src['x'],
                         y = src['y'],
                         psf = src['psf'],
                         amplitude = 1.0,
                         subpix_map = subpix_map)

        src_list.write('%10d' % src_id)
        for var in src_list_vars : src_list.write(' %25.16e' % src[var])
        src_list.write('\n')

    src_list.close()

    hdu_list = pyfits.HDUList([pyfits.PrimaryHDU(image)])
    hdu_list.writeto(image_fname, clobber = True)


