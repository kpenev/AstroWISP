***************
PSF/PRF Fitting
***************

Fitting for Source Shapes
=========================

Fitting for the shapes of point sources (either their PSF or PRF) and an overall
scaling constant is one of the methods of extracting photometry supported by
SuperPhot. This is accomplished using the :mod:`fit_star_shape` module.
Currently only piecewise bi-cubic PSF/PRF models are supported, with the shape
constrained to depend smoothly on image position and any other user-defined
parameters, possibly accross multiple images simultaneously and the amplitudes
(fluxs) of sources being independent of each other. Fitting is done by
constructing an instance of :class:`FitStarShape` and calling it on a collection
of frames to be fit simultaneously and a list of all the sources in each frame.
For details on how to specify fitting parameters and source and frame listts,
see the documentation of :class:`FitStarShape`.

PSF Map utilities
=================
In addition to fitting for the PSF/PRF map and fluxes, 
