"""Collection of utilities for working with files."""

import os
import os.path
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
import logging

from astropy.io import fits

_logger = logging.getLogger(__name__)


def prepare_file_output(
    fname, allow_existing, allow_dir_creation, delete_existing=False
):
    """Ger ready to create/overwrite a file with the given name."""

    result = False
    if os.path.exists(fname):
        if not allow_existing:
            raise OSError(
                f"Destination file {fname!r} already exists and overwritting "
                "not allowed!"
            )
        if delete_existing:
            _logger.info("Overwriting %s", fname)
            os.remove(fname)
        else:
            result = True

    out_path = os.path.dirname(fname)
    if allow_dir_creation and out_path and not os.path.exists(out_path):
        _logger.info("Creating output directory: %s", repr(out_path))
        os.makedirs(out_path)

    return result


@contextmanager
def get_unpacked_fits(fits_fname):
    """Ensure the result is an unpacked version of the frame."""

    with fits.open(fits_fname, "readonly") as fits_file:
        # False positive
        # pylint: disable=no-member
        packed = fits_file[0].header["NAXIS"] == 0
        # pylint: enable=no-member

    if packed:
        with NamedTemporaryFile(
            buffering=0,
            dir=("/dev/shm" if os.path.exists("/dev/shm") else None),
        ) as unpacked_frame:
            with fits.open(fits_fname, "readonly") as fits_file:
                hdu_list = fits.HDUList(
                    [
                        (fits.PrimaryHDU if hdu_ind == 0 else fits.ImageHDU)(
                            header=hdu.header, data=hdu.data
                        )
                        for hdu_ind, hdu in enumerate(fits_file[1:])
                    ]
                )
            hdu_list.writeto(unpacked_frame.name, overwrite=True)
            yield unpacked_frame.name
    else:
        yield fits_fname


def get_fits_fname_root(fits_fname):
    """Return the FITS filename withou directories or extension."""

    result = os.path.basename(fits_fname)
    while True:
        result, extension = os.path.splitext(result)
        if not extension:
            return result


def get_fname_pattern_substitutions(fits_fname, fits_header=None):
    """Return a dictionary that can be used to complete a filename pattern."""

    if fits_header is None:
        with fits.open(fits_fname, "readonly") as fits_file:
            # False positive
            # pylint: disable=no-member
            fits_header = fits_file[
                0 if fits_file[0].header["NAXIS"] else 1
            ].header
            # pylint: enable=no-member

    return dict(
        fits_header,
        FITS_ROOT=get_fits_fname_root(fits_fname),
        FITS_DIR=os.path.dirname(fits_fname),
    )
