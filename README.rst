AstroWISP: ``Astro``\ nomical ``W``\ idefield ``I``\ mages ``S``\ tellar ``P``\ hotometry
=========================================================================================

A tool for extracting stellar photometry from widefild color or monochrome 
images of the night sky

[Full documentation](https://kpenev.github.io/AstroWISP/)

Currently the following photometry methods are supported, and in the future we
plan to add Image Subtraction:

PSF and PRF fitting
-------------------

Fit for the distribution of light from a star on the detector (usually
constrained to vary smoothly across stars) and an individual amplitude for each
star giving a measure of the flux.

Point Sread Function, or PSF, refers to the distribution of light hitting the
detector as a function of the offset from the source center. In order to predict
the value that a particular pixel should have given the PSF one needs to
integrate over the pixel the product of the PSF and the sensitivity of the pixel
at each position within the pixel.

Pixel Response Function, or PRF, incorporates the effect of the detector. The
value of the PRF at a given offset from the source center gives the value a
pixel centered at that location should have. 

AstroWISP allows for both PSF and PRF fitting, imposing a requirement that
either function depends smoothly on the properties of the star being fit. The
dependence is parametrized as an arbitrary polynomial of fnuctions of the
relevant paramaters. Which properties it is allowed to depend on is entirely up
to the user. Typically at least the position of the source center on the
detector is included, but other properties can be included as well (e.g. the
color of the source, temperature of the telescope tube, etc.). Multiple images
can be fit simultaneously imposing the smooth dependence both within and across
images.

Aperture photometry
-------------------

Sum-up the flux in a circular aperture centered around each source. AstroWISP
handles pixels that span the aperture boundary by properly integrating the
product of th PSF and the sub-pixel sensitivity.
