Format for HDF5 Photometry and PSF Fitting Files {#hdf5fileformat_main_page}
================================================

  \brief Describes the format of the HDF5 files produced by FitPSF and
  SubPixPhot.

  Since the most common case for extending photometry files is to add more
  columns as postprocessing steps are performed and it is not necessary to
  add new datapoints non-extendable datasets are used, one for each output
  quantity specified by the user. All files will contain a `Sources` group
  under which all per-source information will be stored. In addition, FitPSF
  output files with global PSF fit will contain an extra group called
  `PSFMap` which will contain the coefficients giving the spatial
  dependence of the PSF on image position. Finally, the `Sources` group
  will have as attributes all the configuration of how it was derived, in
  sufficient detail to exactly reproduce the result.

  The `Sources` group
  --------------------

  ## Attributes: ##

  Most of these are controlled by command line options.

  ### Present for both FitPSF and SubPixPhot: ###

  | Attribute Name | Corresponding Command Line Option/Description |
  | -------------- | --------------------------------------------- |
  | Image          | The filename of the input fits image.         |
  | SubPixMap      | --subpix                                      |
  | PSFModel       | --psf-model                                   |
  | Gain           | --gain                                        |
  | BgAnnulus      | --bg-annulus                                  |
  | InputColumns   | --input-columns                               |
  | FluxMagnitude  | --flux-magnitude                              |
  | IDFormat       | --hat-ids                                     |

  ### SubPixPhot only: ###

  | Attribute Name           | Corresponding Command Line Option                                |
  | ------------------------ | ---------------------------------------------------------------- |
  | AsciiSDKPoly             | --sdk-poly (elliptical Gaussian PSF models only)                 |
  | ExtraError               | --const-error                                                    |
  | Apertures                | --aperture                                                       |
  | MaxExpCoef               | --max-exp-coef                                                   |

  ## Datasets: ##

  Any dataset listed above can be enabled or disabled through the use of
  the `--output-columns` option of the tool producing the file.

  ### Always Available: ###

  | Dataset Name   | Dataset Data Type                           | Output Column Name | 
  | -------------- | ------------------------------------------- | ------------------ |
  | ID             | `string` or `H5T_STD_I16LE`,`H5T_STD_I32LE` | id                 |
  | x              | `H5T_IEEE_F32LE`                            | x                  |
  | y              | `H5T_IEEE_F32LE`                            | y                  |
  | Flux           | `H5T_IEEE_F32LE`                            | flux               |
  | FluxErr        | `H5T_IEEE_F32LE`                            | flux_err           |
  | Magnitude      | `H5T_IEEE_F32LE`                            | mag                |
  | MagnitudeErr   | `H5T_IEEE_F32LE`                            | mag_err            |
  | QualityFlag    | `H5T_C_S1`                                  | flag               |
  | Background     | `H5T_IEEE_F32LE`                            | bg                 |
  | BackgroundErr  | `H5T_IEEE_F32LE`                            | bg_err             |
  | BackgroundNPix | Smallest little endian signed int that fits | nbgpix             |

  ### FitPSF Only: ###

  | Dataset Name     | Dataset Data Type | Output Column Name | 
  | ---------------- | ----------------- | ------------------ |
  | S (only SDK psf) | `H5T_IEEE_F32LE`  | S                  |
  | D (only SDK psf) | `H5T_IEEE_F32LE`  | D                  |
  | K (only SDK psf) | `H5T_IEEE_F32LE`  | K                  |
  | Amplitude        | `H5T_IEEE_F32LE`  | A or amp           |
  | SignalToNoise    | `H5T_IEEE_F32LE`  | sn                 |
  | NFitPixels       | `H5T_STD_I16LE`   | npix               |
  | ReducedChi2      | `H5T_IEEE_F32LE`  | chi2               |

  The `PSFMap` group
  -------------------

  ## Attributes: ##

  | Attribute Name           | Corresponding Command Line Option or Description                                |
  | ------------------------ | ---------------------------------------------------------------- |
  | PSFModel                 | --psf-model                                                      |
  | SubPixMap                | --subpix                                                         |
  | FitOrder                 | --fit-order                                                      |
  | MinS                     | --minS (elliptical Gaussians PSF models only)                    |
  | MaxS                     | --maxS (elliptical Gaussians PSF models only)                    |
  | GridXBoundaries          | --grid (piecewise bicubic PSF models only)                       |
  | GridYBoundaries          | --grid (piecewise bicubic PSF models only)                       |
  | FitTolerance             | --fit-tolerance (elliptical Gaussians PSF models only)           |
  | MaxChi2                  | --max-chi2                                                       |
  | BgExcessThreshold        | --alpha                                                          |
  | MaxSaturatedFraction     | --max-sat-frac                                                   |
  | MinPix                   | --min-pix                                                        |
  | MaxPix                   | --max-pix                                                        |
  | MaxSources               | --max-sources                                                    |
  | MinBgPix                 | --min-bg-pix                                                     |
  | SourceAssignment         | --source-assignment                                              |
  | InitialGuess             | --initial-guess                                                  |
  | MaxAperture              | --max-aperture                                                   |
  | InitialAmplitudeAperture | --initial-amplitude-aperture (piecewise bicubic PSF models only) |
  | MaxAbsAmplitudeChange    | --max-abs-amplitude-change (piecewise bicubic PSF models only)   |
  | MaxRelAmplitudeChange    | --max-rel-amplitude-change (piecewise bicubic PSF models only)   |
  | XOffset                  | Offset of the source x value before evaluating the polynomial    |
  | YOffset                  | Offset of the source y value before evaluating the polynomial    |
  | XScale                   | Scaling of the source x value before evaluating the polynomial   |
  | YScale                   | Scaling of the source y value before evaluating the polynomial   |

  The `PSFMap` group contains only a single dataset (named Coefficients) with
  the datatype being an array of `H5T_IEEE_F32LE`. Each array contains the
  polynomial expansion coefficients of a single PSF parameter. Thus, the
  dataset length is 3 for elliptical Gaussian PSF models and 4 times the
  number of interior grid intersections for piecewise bicubic PSF models.
