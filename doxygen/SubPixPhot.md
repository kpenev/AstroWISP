SubPixPhot {#SubPixPhot_main_page}
==========
  \brief Sub-pixel aperture photometry.

  Implements apperture photometry that uses external knowledge of the PSFs of
  a list of sources and a map of how the sensitivity varies over a single
  detector pixel (all pixels are assumed to have the same sensitivity map) in
  order to accurately account for two effects:
  - pixels only partially intersected by the aperture
  - non-uniform sensitivity of a pixel
  
  The present implementation assumes that the %PSF of each source can be
  described by an elliptical Gaussian on top of a uniform background: \f$ B +
  A\exp\left\{ -\frac{1}{2} \left[S(x^2 + y^2) + D(x^2-y^2) +
  2Kxy\right]\right\}\f$, where x and y are relative to the source center.
  see: <a href="http://arxiv.org/abs/0906.3486"> Pal, A.\ 2009, Ph.D.~Thesis
  </a> for more details
  
  Only the S, D and K parameters (along with the location of the source
  center) need to be supplied on input, since the background and the
  amplitude (or rather the overall flux) of the source are what SubPixPhot
  measures.

  Details of how the corrections are implemented are given
  <a href="SubPixPhot.pdf" target="_blank"><b>here</b></a>

  Different PSF prescriptions can be introduced in the code without actually
  knowing the details of how they are used by the rest of the code. What is
  required is to provide member functions that either calculate integrals
  over rectangles and special pieces of circles, or provide Taylor expansion
  coefficients up to some order. In the latter case, the errors in the
  integrals are empirically estimated (unlike the native case of
  EllipticalGaussianPSF).

  \defgroup SubPixPhot SubPixPhot
  \brief Does aperture photometry correcting for non-uniform pixel
  sensitivity.

  Uses knowledge about the source locations and PSFs on an image, as well as
  an externally provided map for the sensitivity of individual pixels
  (assumed all identical) to calculate the flux within a circular aperture
  around each source. 

  Properly accounts for pixels that stradle the aperture boundary and for the
  fact that some parts of the pixels are more sensitive than others.
