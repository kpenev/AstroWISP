"""Define the :class:`FitStarShape` class, which performs PSF/PRF fitting."""

class FitStarShape:
    """
    Fit for the PSF/PRF of stars and their flux.

    .. _attributes:

    Attributes:
        subpixmap (2D numpy array):    The sub-pixel map, for PSF fitting only.

        shape_terms (str):    The terms the PSF is allowed to depend on. The
            EBNF grammar defining the language for this parameter is::

                (* items in angle brackets (< or >) are assumed to be obvious *)
                (* and thus are not defined. *)

                termchar = <ascii character> - "," - "}" ;

                term = termchar , { termchar } ; (* mathematical expressions  *)
                                                 (* involving variables,      *)
                                                 (* floating point numbers    *)
                                                 (* and pi. The complete list *)
                                                 (* of mathematical functions *)
                                                 (* from c++99's cmath        *)
                                                 (* library are supported.    *)

                list = "{" , term , { "," , term } , "}" ; (* simple listing  *)
                                                           (* of terms to     *)
                                                           (* include.        *)

                poly = "O" , <integer> , list ; (* Expands to all polynomial  *)
                                                (* terms of up to combined    *)
                                                (* order <integer> of the     *)
                                                (* entries in list            *)

                set = list | poly ;

                cross = set , { "*" , set } ; (* expands to the cross product *)
                                              (* of all sets.                 *)

                expression = cross , { \"+\" , cross } ; (* merge the terms   *)
                                                         (* of all cross      *)
                                                         (* products together.*)

        max_chi2 (float):    The value of the reduced chi squared above which
            sources are excluded from the fit. This can indicate non-point
            sources or sources for which the location is wrong among ohter
            things.

        grid (list of floats):    A comma separated list of grid boundaries. Can
            either be a single list, in which case it is used for both the
            horizontal and vertical boundaries. If different splitting is
            desired in the two directions, two lists should be supplied
            separated by ``;``. The first list should contain the vertical (x)
            boundaries and the second list gives the horizontal (y) ones.

        pixel_rejection_threshold (float):    A number defining individual
            pixels to exclude from the PSF fit. Pixels with fitting residuals
            (normalized by the standard deviation) bigger than this value are
            excluded. If zero, no pixels are rejected.

        initial_aperture (float):    This aperture is used to derive an initial
            guess for the amplitudes of sources when fitting for a piecewise
            bicubic PSF model by doing aperture photometry assuming a perfectly
            flat PSF.

        max_abs_amplitude_change (float):    The absolute root of sum squares
            tolerance of the source amplitude changes in order to declare the
            piecewise bicubic PSF fitting converged.

        max_rel_amplitude_change (float):    The relative root of sum squares
            tolerance of the source amplitude changes in order to declare the
            piecewise bicubic PSF fitting converged.

        smoothing (float):    How much smoothing penalty to impose when fitting
            the PSF. ``None`` for no smoothing. Value can be both positive and
            negative and will always result in smoothing (less for negative
            values).

        min_convergence_rate (float):    If the rate of convergence falls below
            this threshold, iterations are stopped. The rate is calculated as
            the fractional decrease in the difference between the amplitude
            change and the value when it would stop, as determined by the
            :attr:`max_abs_amplitude_change` and
            :attr:`max_rel_amplitude_change` attributes.

        max_iterations (int):    No more than this number if iterations will be
            performed. If convergence is not achieved before then, the latest
            estimates are output and an exception is thrown. A negative value
            allows infinite iterations. A value of zero, along with an initial
            guess for the PSF causes only the amplitudes to be fit for PSF
            fitting photometry with a known PSF. It is an error to pass a value
            of zero for this option and not specify and initial guess for the
            PSF.
    """

    def __init__(self,
                 *,
                 subpixmap,
                 shape_terms,
                 max_chi2,
                 grid,
                 pixel_rejection_threshold,
                 initial_aperture,
                 max_abs_amplitude_change,
                 max_rel_amplitude_change,
                 smoothing,
                 min_convergence_rate,
                 max_iterations):
        """
        Set-up an object ready to perform PSF/PRF fitting.

        Args:
            keyword_only:    See class attributes_.

        Returns:
            None
        """

        self.subpixmap = subpixmap
        self.shape_terms = shape_terms
        self.max_chi2 = max_chi2
        self.grid = grid
        self.pixel_rejection_threshold = pixel_rejection_threshold
        self.initial_aperture = initial_aperture
        self.max_abs_amplitude_change = max_abs_amplitude_change
        self.max_rel_amplitude_change = max_rel_amplitude_change
        self.smoothing = smoothing
        self.min_convergence_rate = min_convergence_rate
        self.max_iterations = max_iterations

    def __call__(self, image_sources):
        """
        Fit for the shape of the sources in a collection of imeges.

        Args:
            image_sources (dict):    Dictionary, with keys - the filename of the
                reduced frame (FITS) to fit the PSF/PRF of and values - list of
                sources to process, defining at least the following quantities:

                    * **ID** (string): some unique identifier for the source

                    * **x** (float): The x coordinate of the source center in pixels

                    * **y** (float): See ``x``

                May define additional quantities on which the PSF shape is
                allowed to depend. This can be either a numy record array with
                field names as keys or a dictionary with field names as keys
                and 1-D numpy arrays of identical lengths as values.

        Returns:
            dict:

                * keys:
                    The filenames of the reduced frames used in fitting.

                * values:
                    2-tuples containing the following arrays:

                    * 4D numpy array:
                        The coefficients of the PSF/PRF map.

                    * 1D numpy array:
                        The best-fit fluxes of the sources in the same order as
                        specified in the :obj:`image_sources` argument.
        """

        raise Exception('Not implemented')
