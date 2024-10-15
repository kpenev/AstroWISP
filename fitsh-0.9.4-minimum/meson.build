#Control file for meson/meson-python (compiling and PIP packaging)

project(
  'fitsh-astrowisp',
  'c',
  version: '0.0.0',
  meson_version: '>= 1.4.0',
  default_options: [
    'optimization=3',
    'c_std=c11',
    'prefer_static=true'
  ]
)

dependencies = []

add_project_arguments(
  '-D_GNU_SOURCE',
  '-D_FILE_OFFSET_BITS=64',
  '-D_LARGEFILE64_SOURCE',
  '-DSCANARG_WITHOUT_FNMATCH',
  '-D_FITSH_SOURCE',
  language: 'c'
  )

if host_machine.system() == 'windows'
  add_project_arguments(
    '-D_USE_MATH_DEFINES',
    '-DNO_VISIBILITY',
    '-static',
    '-static-libgcc',
    language: 'c'
    )
else
  add_project_arguments(
    '-fvisibility=hidden',
    '-DNIX_VISIBILITY',
    language : 'c'
  )
endif

includedir = include_directories('../include')

subdir('libfits')
subdir('src')
