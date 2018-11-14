"""A collection of functions used by the fit_star_shape unit tests."""

from math import ceil
import os.path
import sys
from ctypes import c_double
import numpy

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '..',
            '..'
        )
    )
)

#Needs to be after sys.path is updated to allow importing superphot package.
#pylint: disable=wrong-import-position
from superphot.fake_image.image import Image
#pylint: enable=wrong-import-position

def make_image_and_source_list(sources,
                               extra_variables,
                               subpix_map):
    """
    Create an image and a list of the sources in it ready for psf fitting.

    Args:
        sources:    A list of dictionaries with at least the following keywords:

                * x:    The x coordinate of the source center.

                * y:    The y coordinate of the source center.

                * psf:    An instance of some sub-class of PSFBase giving the
                  sources's PSF. It should already be scaled to the desired
                  flux.

            Additional keywords may be added to the source list and hence
            available as variables for PSF fitting by listing the names in the
            extra_variables argument.

        extra_variables:    A list of additional keywords from sources to add to
            the source list and the order in which those should be added. The
            corresponding entries in sources must be floating point values.

        subpix_map:    The sub-pixel map to impose on the image. For more
            details see same name argument of Image.add_source.

    Returns:
        numpy record array:
            The sources added to the image. The fields give the variables
            defined for the sources.
    """

    min_x = min(s['x'] for s in sources)
    max_x = max(s['x'] for s in sources)
    min_y = min(s['y'] for s in sources)
    max_y = max(s['y'] for s in sources)

    image = Image(int(ceil(max_x + min_x)),
                  int(ceil(max_y + min_y)),
                  background=1.0)

    src_list_vars = ['x', 'y'] + extra_variables
    src_list = numpy.empty(
        len(sources),
        dtype=(
            [('ID', '|S6')]
            +
            [(var, c_double) for var in src_list_vars]
        )
    )


    for src_id, src in enumerate(sources):
        image.add_source(x=src['x'],
                         y=src['y'],
                         psf=src['psf'],
                         amplitude=1.0,
                         subpix_map=subpix_map)

        src_list[src_id]['ID'] = b'%06d' % src_id
        for var in src_list_vars:
            src_list[src_id][var] = src[var]

    return image, src_list
