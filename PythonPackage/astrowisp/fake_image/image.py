"""Defines a base class for fake star images."""

from math import floor, ceil
import numpy

#This inherits a whole lot of methods from numpy.ndarray, adding 2 extra
#justifies the inherited class.
#pylint: disable=too-few-public-methods
class Image(numpy.ndarray):
    """A base class for fake images."""

    def __new__(cls, x_resolution, y_resolution, background=0.0):
        """Create an image with the given resolution and background."""

        result = super().__new__(cls,
                                 shape=(y_resolution, x_resolution),
                                 dtype=float)
        result.fill(background)
        return result

    #x and y are perfectly reasonable names for 2-D position coordinates.
    #pylint: disable=invalid-name
    #
    #Not counting self and only allowing psf or prf to be specified, there are
    #really only 5 arguments.
    #pylint: disable=too-many-arguments
    #
    #Removing local variables makes the function less readable
    #pylint: disable=too-many-locals
    def add_source(self,
                   x,
                   y,
                   amplitude,
                   psf=None,
                   prf=None,
                   subpix_map=numpy.ones((1, 1))):
        """
        Adds a source to this image.

        Args:
            x:    The x coordinate of the center of the new source relative to
                the lower left corner of the image.

            y:    The x coordinate of the center of the new source relative to
                the lower left corner of the image.

            amplitude:    The amplitude to scale the PSF by.

            psf:    The PSF to use for the new source. Should be of some type
                inherited from fake_image.psf or None if PRF is used instead.

            prf:    The pixel response function (PRF) to use for the new
                source. Should be of some type inherited from fake_image.psf or
                None if PSF is used instead. The pixel response function is the
                PSF convolved with the pixel sensitivity map.

            subpix_map:    The sub-pixel sensitivity map. Should be a
                numpy.ndarray of some sort (i.e. provide shape attribute and
                subscripting).

        Returns:
            None
        """

        y_res, x_res = self.shape
        for pix_y in range(max(0, floor(y - psf.get_down_range())),
                           min(y_res, ceil(y + psf.get_up_range()))):
            for pix_x in range(max(0, floor(x - psf.get_left_range())),
                               min(x_res, ceil(x + psf.get_right_range()))):
                if prf is None:
                    assert subpix_map is not None
                    y_splits, x_splits = subpix_map.shape
                    for subpix_y in range(y_splits):
                        for subpix_x in range(x_splits):
                            self[pix_y, pix_x] += (
                                amplitude
                                *
                                subpix_map[subpix_y, subpix_x]
                                *
                                psf.integrate(pix_x + subpix_x / x_splits - x,
                                              pix_y + subpix_y / y_splits - y,
                                              1.0 / x_splits,
                                              1.0 / y_splits)
                            )
                else:
                    assert psf is None
                    self[pix_y, pix_x] += amplitude * prf(pix_x - x + 0.5,
                                                          pix_y - y + 0.5)
    #pylint: enable=invalid-name
    #pylint: enable=too-many-arguments
    #pylint: enable=too-many-locals
#pylint: enable=too-few-public-methods
