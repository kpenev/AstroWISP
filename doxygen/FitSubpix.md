FitSubpix   {#FitSubpix_main_page}
=========
  \brief Derives a best-fit sub-pixel sensitivity map.

  Fits for the sub-pixel sensitivity of pixels that minimizes the variance of
  the flux measured for a set of sources from a number of input images.
  Information on the location of the sources as well as their %PSF on each
  image is required. Several fitting methods are implemented:
  - Simplex: uses the GNU Scientific Library (GSL) simplex minimization.
  - Newton-Raphson: uses analytical expressions for the first and second
	derivatives of the flux derived for each source with respect to the
	individual sub-pixel sensitivity values. Details are given
	<a href="NRFitting.pdf" target="_blank"><b>here</b></a>
  - MultiNest: uses the MultiNest algorithm (as implemented in MultiNest
	v2.14).
  - Simulated annealing: uses the GSL simulated-annealing implementation with
	parameters that can be tuned on the command line.

  The ideal input dataset consists of many bright but not saturated sources
  with a %PSF that changes significantly on the scale of one pixel, the
  locations of which move around from image to image by one or several
  pixels.

\defgroup FitSubpix FitSubpix
