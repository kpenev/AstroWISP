py = import('python').find_installation(pure: false)

psf_lib = static_library(
  'psf',
  ['psf-io.c', 'psf-base.c', 'psf-determine.c'],
  include_directories: includedir,
  install: false
)

common_src = [
  'statistics.c',
  'ui.c',
  'longhelp.c',
  'io/scanarg.c',
  'io/iof.c',
  'io/tokenize.c'
]

fistar_src = [
  'fistar.c',
  'math/fit/lmfit.c',
  'math/poly.c',
  'math/polyfit.c',
  'math/spline/biquad.c',
  'math/spline/biquad-isc.c',
  'math/expint/expint.c',
  'math/intersec/intersec-cri.c',
  'math/delaunay.c',
  'index/sort.c',
  'index/multiindex.c',
  'link/floodfill.c',
  'link/linkpoint.c',
  'link/linkblock.c',
  'fistar-io.c',
  'fitsmask.c',
  'tensor.c',
  'star-base.c', 
  'star-model.c', 
  'star-psf.c',
  'star-draw.c',
  'star-cand-pp.c',
  'star-cand-biq.c',
  'star-cand-trb.c',
  'star-cand-lnk.c', 
  'common.c',
  'imgtrans.c',
  'background.c',
  'magnitude.c'
]

grcollect_src = [
  'grcollect.c',
  'cache.c',
  'str.c'
]

math_dependence = meson.get_compiler('c').find_library('m', required : false)

executable(
  'fistar',
  fistar_src + common_src,
  include_directories: includedir,
  link_with: [psf_lib, fits_lib],
  dependencies: math_dependence,
  install: true,
  install_dir: py.get_install_dir() / 'astrowisp'
)

executable (
  'grcollect',
  grcollect_src + common_src,
  include_directories: includedir,
  link_with: fits_lib,
  dependencies: math_dependence,
  install: true,
  install_dir: py.get_install_dir() / 'astrowisp'
)
