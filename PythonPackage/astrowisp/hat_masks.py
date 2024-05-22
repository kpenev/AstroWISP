"""Utilities for working with HAT-style masks in FITS headers."""

import numpy

from astrowisp._initialize_library import get_astrowisp_library

mask_flags = dict(get_astrowisp_library().mask_flags)

def parse_hat_mask(header):
    """
    Extract the HAT-style mask contained in header.

    Args:
        header:    The header of the image whose mask to parse.

    Returns:
        numpy.array(dtype=uint8):
            array with exactly the same resolution as the input image containing
            a bit-field for each pixel indicating any bad-pixel flags raised per
            the header.

    Examples:

        >>> from astropy.io import fits
        >>>
        >>> with fits.open('/Users/kpenev/tmp/1-447491_4.fits.fz',
        >>>                mode='readonly') as f:
        >>>     image_mask = parse_hat_mask(f[1].header)
        >>>
        >>>     flag_name = 'OVERSATURATED'
        >>>
        >>>     matched = numpy.bitwise_and(image_mask,
        >>>                                 mask_flags[flag_name]).astype(bool)
        >>>
        >>>     #Print number of pixels for which the OVERSATURATED flag is
        >>>     #raised
        >>>     print(flag_name + ': ' + repr(matched.sum()))
        >>>
        >>>     #Output x, y, flux for the pixels flagged as OVERSATURATED
        >>>     for y, x in zip(*numpy.nonzero(matched)):
        >>>         print('%4d %4d %15d' % (x, y, f[1].data[y, x]))
    """

    mask_string = ''.join((c[1] + ' ') if c[0] == 'MASKINFO' else ''
                          for c in header.items()).encode('ascii')
    mask = numpy.zeros((header['NAXIS2'], header['NAXIS1']), dtype='int8')
    get_astrowisp_library().parse_hat_mask(mask_string,
                                           header['NAXIS1'],
                                           header['NAXIS2'],
                                           mask.ravel())
    return mask
