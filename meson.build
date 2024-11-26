#Control file for meson/meson-python (compiling and PIP packaging)

project(
    'astrowisp',
    'cpp', 'c',
    version: '1.0.1',
    meson_version: '>= 1.4.0',
    default_options: [
        'optimization=3',
        'cpp_std=c++14',
        'prefer_static=true'
    ]
)

py = import('python').find_installation(pure: false)

compiler_id = meson.get_compiler('cpp').get_id()
toolchain = {
    'gcc': 'TOOLCHAIN_GCC',
    'clang': 'TOOLCHAIN_CLANG',
    'msvc': 'TOOLCHAIN_MSVC'
}[compiler_id]

add_project_arguments('-D' + toolchain,
                      language : 'cpp')

add_project_arguments(
  '-D_GNU_SOURCE',
  '-D_FILE_OFFSET_BITS=64',
  '-D_LARGEFILE64_SOURCE',
  '-D_FITSH_SOURCE',
  language: 'c'
  )


if get_option('buildtype') == 'debug'
    add_project_arguments('-DVERBOSE_DEBUG', '-DDEBUG', language: 'cpp')
elif get_option('buildtype') == 'release'
    add_project_arguments('-DNDEBUG', language: 'cpp')
endif

dependencies = [
    dependency('boost', modules: ['program_options']),
    dependency('eigen3')
]


if host_machine.system() == 'windows'
    add_project_arguments(
      '-DBOOST_RESULT_OF_USE_TR1',
      '-D_USE_MATH_DEFINES',
      '-DNO_VISIBILITY',
      '-static',
      '-static-libgcc',
      '-static-libstdc++',
      language: 'cpp'
    )

  add_project_arguments(
    '-D_USE_MATH_DEFINES',
    '-DNO_VISIBILITY',
    '-DHOST_WIN32',
    '-DSCANARG_WITHOUT_FNMATCH',
    '-static',
    '-static-libgcc',
    language: 'c'
  )

else
    add_project_arguments(
      '-fvisibility=hidden',
      '-DNIX_VISIBILITY',
      language : 'cpp'
    )

    add_project_arguments(
      '-fvisibility=hidden',
      '-DNIX_VISIBILITY',
      language : 'c'
    )

endif

cpp_src = [
    'src/Core/ParseCSV.cpp',
    'src/Core/SDKSource.cpp',
    'src/Core/SubPixelMap.cpp',
    'src/Core/CInterface.cpp',
    'src/Core/CommandLineUtil.cpp',
    'src/Core/SourceID.cpp',
    'src/IO/parse_hat_mask.cpp',
    'src/IO/OutputArray.cpp',
    'src/IO/H5IODataTree.cpp',
    'src/IO/Binary.cpp',
    'src/IO/CInterface.cpp',
    'src/IO/parse_grid.cpp',
    'src/IO/CommandLineUtil.cpp',
    'src/IO/CommandLineConfig.cpp',
    'src/Background/Manual.cpp',
    'src/Background/Annulus.cpp',
    'src/Background/Source.cpp',
    'src/Background/CInterface.cpp',
    'src/Background/MeasureAnnulus.cpp',
    'src/Background/CommandLineUtil.cpp',
    'src/Background/MeasureExcludingSources.cpp',
    'src/PSF/PSF.cpp',
    'src/PSF/Map.cpp',
    'src/PSF/IntegralUtil.cpp',
    'src/PSF/EllipticalGaussian.cpp',
    'src/PSF/MapSourceContainer.cpp',
    'src/PSF/EllipticalGaussianIntegralWedge.cpp',
    'src/PSF/Piecewise.cpp',
    'src/PSF/PiecewiseBicubicCell.cpp',
    'src/PSF/PiecewiseBicubicMap.cpp',
    'src/PSF/EllipticalGaussianIntegralByOrder.cpp',
    'src/PSF/EllipticalGaussianMap.cpp',
    'src/PSF/CInterface.cpp',
    'src/PSF/EllipticalGaussianIntegralRectangle.cpp',
    'src/PSF/WedgeIntegral.cpp',
    'src/PSF/DataTreeCalculations.cpp',
    'src/PSF/CirclePieceIntegral.cpp',
    'src/PSF/CommandLineUtil.cpp',
    'src/PSF/LocalPolynomial.cpp',
    'src/FitPSF/LinearSource.cpp',
    'src/FitPSF/IOSources.cpp',
    'src/FitPSF/CInterface.cpp',
    'src/FitPSF/PiecewiseBicubicPSFSmoothing.cpp',
    'src/FitPSF/Common.cpp',
    'src/FitPSF/Pixel.cpp',
    'src/FitPSF/Config.cpp',
    'src/FitPSF/PiecewiseBicubic.cpp',
    'src/SubPixPhot/CInterface.cpp',
    'src/SubPixPhot/Config.cpp',
    'src/SubPixPhot/SourceInput.cpp',
    'src/SubPixPhot/SubPixPhotIO.cpp'
]

shared_library(
    'astrowisp',
    cpp_src,
    dependencies: dependencies,
    install: true,
    install_dir: py.get_install_dir() / 'astrowisp'
)

subdir('PythonPackage/astrowisp')
subdir('fitsh-0.9.4-minimum')
