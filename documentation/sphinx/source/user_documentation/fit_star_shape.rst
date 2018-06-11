***************
PSF/PRF Fitting
***************

Fitting for the shapes of point sources (either their PSF or PRF) and an overall
scaling constant is one of the methods of extracting photometry supported by
SuperPhot. This is accomplished using the :mod:`fit_star_shape` module.
Currently only piecewise bi-cubic PSF/PRF models are supported, with the shape
constrained to depend smoothly on image position and any other user-defined
parameters, possibly accross multiple images simultaneously and the amplitudes
(fluxs) of sources being independent of each other. Fitting is done by
constructing an instance of :class:`FitStarShape` and calling it on a collection
of frames to be fit simultaneously and a list of all the sources in each frame.
The exact 
accomplished by
specifying the frames to process and a list of sources for each frame
The following configuration
parameters control the fitting:

    * **Source lists:** dictionary, with keys - the filename of the reduced frame to
      fit the PSF/PRF of and values - list of sources to process, defining at
      least the following quantities:

        * **ID** (string): some unique identifier for the source

        * **x** (float): The x coordinate of the source center in pixels

        * **y** (float): See ``x``

      May define additional quantities on which the PSF shape is allowed to
      depend. This can be either a numy record array with field names as keys or
      a dictionary with field names as keys and 1-D numpy arrays of identical
      lengths as values.

As usual those can be specified in three stages: when constructing the instance
of 
