c_src = [
  'fits-common.c',
  'fits-header.c',
  'fits-image.c',
  'fits-table.c',
  'fits-bintable.c',
  'fits-draw.c',
  'fits-ui.c',
  'fits-core.c'
]

fits_lib = static_library(
  'libfits',
  c_src,
  include_directories: includedir,
  install: false
)
