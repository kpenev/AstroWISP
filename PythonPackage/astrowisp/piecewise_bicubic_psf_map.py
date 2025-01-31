"""A wrapper class for working with PSF/PRF maps from the C/C++ library."""

from astrowisp.piecewise_bicubic_psf import PiecewiseBicubicPSF
from astrowisp._initialize_library import get_astrowisp_library

class PiecewiseBicubicPSFMap:
    """Provide convenient python interface to shape fitting results."""

    def __init__(self, star_shape_map_tree):
        """
        Prepare to query the map generated by a star shape fit.

        Args:
            star_shape_map_tree(IOTree):    The result returned by
                calling FitStarShape.fit().

        Returns:
            None
        """

        self._astrowisp_library = get_astrowisp_library()
        self._library_map = (
            self._astrowisp_library.create_piecewise_bicubic_psf_map(
                star_shape_map_tree.library_tree
            )
        )

    def __call__(self, term_values):
        """
        Evaluate the map for a given set of terms.

        Args:
            term_values:    The terms that PSF parameters depend on evaluated
                for the particular source we wish to know the PSF of.

        Returns:
            PSF:
                The PSF/PRF the map predicts for the given arguments.
        """

        assert (
            (len(term_values.shape) == 2 and term_values.shape[0] == 1)
            or
            len(term_values.shape) == 1
        )
        term_values = term_values.flatten()
        return PiecewiseBicubicPSF(
            self._astrowisp_library.evaluate_piecewise_bicubic_psf_map(
                self._library_map,
                term_values
            )
        )

    def __del__(self):
        """Delete any objects allocated by the library."""

        self._astrowisp_library.destroy_piecewise_bicubic_psf_map(
            self._library_map
        )
