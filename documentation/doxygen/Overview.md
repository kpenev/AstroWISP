Sub-Pixel Photometry	{#mainpage}
====================

This project provides three tools, each addressing a a different task related
to doing aperture photometry in an optimal way for the case typical of ground
based transit searches for extrasolar planets. The emphasis is on correcting
for the fact that the sensitivity of some detectors is not uniform over a
pixel. This effect will cause different flux to be measured for the same
source depending on where within a pixel its center of light falls.

The three tools are:

 - \ref SubPixPhot_main_page - performs aperture photometry on a list of sources with known
   locations and point spread functions (PSFs), correcting for:
   - pixels only partially intersected by the aperture
   - non-uniform sensitivity of a pixel (a sub-pixel sensitivity map must be
     provided)
 - \ref FitSubpix - tries to  fit for the sub-pixel sensitivity from a set of
   images which must contain the same sources with the same brightness (at
   least that's the ideal).
 - \ref FitPSF - fits for the PSFs of a given list of sources on an image either
   on a source-by-source basis, or as an ensamble under the assumption that
   the PSF varies only slowly over the image. If a subpixel sensitivity map
   is provided it is used when fitting for the %PSF.

The sub-pixel sensitivity map is assumed to be constant on rectangular
pieces. It is supplied as a fits image.

