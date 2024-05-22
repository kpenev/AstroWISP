*******************
Aperture Photometry
*******************

A commonly used photometry method is to sum the flux within a circular aperture
centered on each source. Because large apertures include more sky noise which
affects faint sources, and small apertures ignore flux from bright sources, thus
increasing the Poisson noise, multiple apertures are required to get optimal
photometry for all sources. Two sources of error can creep into this
measurement. First, for pixels entirely within the aperture, non--uniform pixel
sensitivity can align differently with the non--uniform PSF from one image to
another producing a different response.  Second, some pixels inevitably cross
the aperture boundary. In addition to the former effect, the fraction of the
flux that should be counted as inside the aperture for those pixels is in
general *not equal to* the fraction of the pixel that is inside the aperture.
The :class:`astrowisp.SubPixPhot` tool handles both of these effects exactly.
For each pixel, the pixel response is multiplied by the integral of the PSF over
the part of the pixel inside the aperture and divided by the integral of the
product of the PSF and the pixel sensitivity function over the entire pixel
before being added to the total flux. Both integrals of this procedure are
calculated analytically and without any approximations.

Aperture photometry with :class:`astrowisp.SubPixPhot` is demonstrated `here
<example_mock_data.ipynb#Aperture-Photometry>`_
