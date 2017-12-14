FitPSF  {#FitPSF_main_page}
======
  \brief Derives best fit PSFs correcting for non-uniform pixel sensitivity.

  Non-uniform sensitivity on a sub-pixel level affects the %PSF in an
  image, and the effect depends on where within a pixel the center of the
  source lies. Thus, it is necessary to correct for the sub-pixel sensitivity
  variations when deriving PSFs for a given image. The FitPSF tool does just
  that.

  Since we envision that %PSF fitting will be performed on images taken with
  a small telescope with pixels that are quite large on the sky, locations of
  the centers of sources can be determined much more accurately by finding an
  astrometric transformation that maps an external high resolution catalogue
  (e.g. 2MASS, GSC, etc.) onto the images than by measuring the center of
  light of individual sources. For this reason, FitPSF requires that an
  external list of source positions is provided, for which fitting the %PSF
  shape and amplitude is necessary.

  %PSF fitting begins by assigning pixels to each input source. The details
  of how that is done depend on the PSF model. At present two PSF models are
  available: <a href="PSFFitting.pdf" target="_blank"> elliptical
  gaussians </a> and
  @ref PiecewiseBicubicFitPSF_main_page "piecewise bicubics". After
  source pixels are identified, the residual between predicted detector
  counts for each pixel and the measured counts is minimized. 

\defgroup FitPSF FitPSF
