project('superphot', 
        'cpp',
        default_options: ['optimization=3', 'cpp_std=c++14'])

src_files = ['src/Core/ParseCSV.cpp', 
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
             'src/SubPixPhot/SubPixPhotIO.cpp']
compiler_id = meson.get_compiler('cpp').get_id()
toolchain = {
  'gcc': 'TOOLCHAIN_GCC', 
  'clang': 'TOOLCHAIN_CLANG',
  'msvc': 'TOOLCHAIN_MSVC'
}[compiler_id]
dependencies = [
  dependency('boost', modules: ['program_options']),
  dependency('eigen3')
]
add_project_arguments('-D' + toolchain, language : 'cpp')
if compiler_id == 'msvc'
  add_project_arguments('-D_USE_MATH_DEFINES', language: 'cpp')
  add_project_arguments('-DBOOST_RESULT_OF_USE_TR1', language: 'cpp')
elif 
  add_project_arguments('-fvisibility=hidden', language : 'cpp')
endif


library('superphot', src_files, dependencies: dependencies)
