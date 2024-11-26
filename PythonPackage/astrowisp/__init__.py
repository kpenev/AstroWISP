"""AstroWISP python interface."""

from os import path

from astrowisp.background import BackgroundExtractor
from astrowisp.fit_star_shape import FitStarShape
from astrowisp.subpixphot import SubPixPhot
from astrowisp.io_tree import IOTree
from astrowisp.piecewise_bicubic_psf_map import PiecewiseBicubicPSFMap
from astrowisp.piecewise_bicubic_psf import PiecewiseBicubicPSF

_module_path = path.dirname(path.abspath(__file__))

fistar_path = path.join(_module_path, 'fistar')
if not path.exists(fistar_path):
    fistar_path += '.exe'
assert path.exists(fistar_path)

__all__ = ['BackgroundExtractor',
           'FitStarShape',
           'SubPixPhot',
           'PiecewiseBicubicPSF',
           'PiecewiseBicubicPSFMap']
