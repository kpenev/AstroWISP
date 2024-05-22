"""A collection of functions used by the fit_star_shape unit tests."""

from math import ceil
from ctypes import c_double, c_bool
import numpy

from asteval import Interpreter

from astrowisp.fake_image.image import Image

def make_image_and_source_list(sources,
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

                * enabled(optional):    True/False flag of whether to include
                    this source in shape fitting.

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

    src_list_dtype = (
        [('ID', '|S6')]
        +
        [(var, c_double) for var in 'xy']
    )
    if 'enabled' in sources[0]:
        src_list_dtype.append(('enabled', c_bool))

    src_list = numpy.empty(
        len(sources),
        dtype=src_list_dtype
    )


    for src_id, src in enumerate(sources):
        image.add_source(x=src['x'],
                         y=src['y'],
                         psf=src['psf'],
                         amplitude=1.0,
                         subpix_map=subpix_map)

        src_list[src_id]['ID'] = b'%06d' % src_id
        for var in src_list.dtype.names[1:]:
            src_list[src_id][var] = src[var]

    return image, src_list

def evaluate_psffit_terms(sources, terms_str):
    """
    Evaluate the specified terms to properly formatted array for PSF fitting.

    Args:
        sources:    See `sources` argument to `make_image_and_source_list()`.

        terms_str([str]):    List of strings to evaluate using entries in
            sources.

    Returns:
        2-D numpy array:
            The terms the PSF map is allowed to depend on evaluated for each
            source, organized suitably to pass as the fifth argument to
            `FitStarShape.fit()`.
    """

    num_sources = len(sources)
    num_terms = len(terms_str)

    result = numpy.empty((num_sources, num_terms), dtype=c_double)

    evaluate = Interpreter()
    for src_ind, src in enumerate(sources):
        evaluate.symtable.update(src)
        for term_ind, term in enumerate(terms_str):
            result[src_ind, term_ind] = evaluate(term)

    #print('PSF fit terms: ' + repr(result))

    return result
