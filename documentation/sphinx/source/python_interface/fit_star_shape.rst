***************
PSF/PRF Fitting
***************

Fitting for Source Shapes
=========================

Fitting for the shapes of point sources (either their PSF or PRF) and an overall
scaling constant is one of the methods of extracting photometry supported by
AstroWISP. This is accomplished using the :mod:`astrowisp.fit_star_shape`
module.  Currently only piecewise bi-cubic PSF/PRF models are supported, with
the shape constrained to depend smoothly on image position and any other
user-defined parameters, possibly accross multiple images simultaneously and the
amplitudes (fluxes) of sources being independent of each other. Fitting is done
by constructing an instance of :class:`astrowisp.FitStarShape` and calling it on
a collection of frames to be fit simultaneously and a list of all the sources in
each frame.  For details on how to specify fitting parameters and source and
frame listts, see the documentation of :class:`astrowisp.FitStarShape`.

An example of fitting for the PSF in an image can be found
`here <example_mock_data.ipynb#PSF-fitting>`_

PSF Map utilities
=================

Smooth dependence of PSF/PRF on parameters is enforced by modeling the PSF/PRF
parameters as low order polynomials of user specified quantities. AstroWISP
defines a language for specifying the dependence (see
:attr:`astrowisp.FitStarShape.shape_terms` for the definition of the language).
The :class:`astrowisp.SmoothDependence` class offers tools for parsing the
language, and generating the various terms involved in the dependence.
