********************
Measuring Background
********************

Photometry starts by determining the background level and its uncertainty for
each source. At present AstroWISP uses the median average and the scatter around
it (after outlier rejection) of pixels in an annulus around each source,
excluding pixels too close to other sources. The user specifies two radii:

    * Inner radius, defining pixels that are so close to a source that they
      include not just background, but source light. Pixels with centers within
      this radius around each source do not contribute to any source's
      background value.

    * Outer radius, defining pixels local enough to the source to represent the
      background under it. Each source's background is then the median of all
      pixels within this outer radius that are not within the the inner radius
      of this source or any other source.

Background measurement is performed using the
:class:`astrowisp.BackgroundExtractor` class and requires specifying the image,
list of source locations, the inner/outer background radii  and returns
estimates of the background for each source, its error (standard deviation), and
the number of pixels that participated in the determination of the background
value and error::

    from astrowisp import BackgroundExtractor
    import numpy.random

    image = numpy.random.rand(100, 100)
    x = numpy.random.rand(10) * 1000.0
    y = numpy.random.rand(10) * 1000.0

    measure_background = BackgroundExtractor(image,
                                             inner_radius=6.0,
                                             outer_radius=13.0)
    bg_value, bg_error, bg_num_pix = measure_background(x, y)
