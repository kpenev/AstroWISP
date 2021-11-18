/**\file
 *
 * \brief The implementation of the FitPSF::Config methods.
 *
 * \ingroup FitPSF
 */

#include "Config.h"

namespace FitPSF {

    void Config::describe_options()
    {

        opt::options_description io_options("Input/Output options");
        io_options.add_options()
            (
                "io.source-list,l",
                opt::value<std::string>()->default_value(""),
                "File containing the input source list."
            )
            (
                "io.subpix,s",
                opt::value<std::string>()->default_value(""),
                "FITS file representing the sensitivity structure of a "
                "pixel. Uniform sensitivity is assumed if empty string. This"
                " option is ignored if PRF fitting is performed."
            )
            (
                "io.initial-guess",
                opt::value<std::string>()->default_value(""),
                "A file to read an initial guess for the PSF map from."
            )
            (
                "io.input-columns",
                opt::value< Core::StringList >(),
                "A comma separated list of input column names. Recognized "
                "values are the same as in output-quantities."
            )
            (
                "io.expect-error-hdu",
                opt::value<unsigned>(),
                "If this option is non-zero, estimated errors in the pixel "
                "values of the primary image are assumed to be contained in the"
                " HDU with the given number."
            )
            (
                "io.error-image",
                opt::value<std::string>()->default_value(""),
                "Use the primary HDU in the given FITS file as the error "
                "estimate for the pixel fluxes. This option supercedes "
                "io.expect-error-hdu."
            );

        opt::options_description psffit_options("Options defining the fit");
        psffit_options.add_options()
            (
                "psf.model,m",
                opt::value<PSF::ModelType>(),
                "The PSF/PRF model to fit. Case insensitive. One of: sdk, "
                "bicubic."
            )
            (
                "psf.max-chi2",
                opt::value<double>(),
                "The value of the reduced chi squared above which sources "
                "are excluded from the fit. This can indicate non-point "
                "sources or sources for which the location is wrong among "
                "ohter things. "
            )
            (
                "psf.min-convergence-rate",
                opt::value<double>(),
                "If the rate of convergence falls below this threshold, "
                "iterations are stopped. For piecewise bicubic PSF/PRF "
                "fitting, the rate is calculated as the fractional decrease "
                "in the difference between the amplitude change and the "
                "value when it would stop, as determined by the "
                "--psf.bicubic.max-abs-amplitude-change and "
                "--psf.bicubic.max-rel-amplitude-change parameters. For SDK "
                "PSF/PRF fitting <++>."
            )
            (
                "psf.max-iterations",
                opt::value<int>(),
                "No more than this number if iterations will be performed. "
                "If convergence is not achieved before then, the latest "
                "estimates are output and an exception is thrown. A negative"
                " value allows infinite iterations. A value of zero, along "
                "with an initial guess for the PSF/PRF causes only the "
                "amplitudes to be fit for PSF/PRF fitting photometry with a "
                "known PSF/PRF. It is an error to pass a value of zero for "
                "this option and not specify and initial guess for the "
                "PSF/PRF."
            )
            (
                "psf.ignore-dropped",
                opt::value<bool>(),
                "If this option is passed, sources dropped during source "
                "selection will not have their amplitudes fit for. Instead "
                "their PSF fit fluxes/magnitudes and associated errors will "
                "all be NaN. This is useful if PSF fitting is done solely "
                "for the purpose of getting a PSF map."
            )
            (
                "psf.sdk.minS",
                opt::value<double>(),
                "The minimum value for the S PSF/PRF paramater to consider "
                "for the fit. Ignored unless the PSF model is sdk."
            )
            (
                "psf.sdk.maxS",
                opt::value<double>(),
                "The maximum value for the S PSF/PRF paramater to consider "
                "for the fit. You need to be careful with this. Severe "
                "overestimation will result in painfully slow PSF fitting. "
                "Ignored unless the PSF model is sdk."
            )
            (
                "psf.sdk.fit-tolerance",
                opt::value<double>(),
                "The required size of the simplex before a fit is declared "
                "converged. Ignored for bicubic PSF/PRF model."
            )
            (
                "psf.sdk.use-simplex",
                opt::value<bool>(),
                "If this options is passed, fitting will be done using the "
                "GSL simplex minimizer instead of the much faster "
                "Newton-Raphson. Ignored for the bicubic PSF/PRF model."
            )
            (
                "psf.bicubic.grid",
                opt::value<PSF::Grid>(),
                "A comma separated list of grid boundaries. Can either be a "
                "single list, in which case it is used for both the "
                "horizontal and vertical boundaries. If different splitting "
                "is desired in the two directions, two lists should be "
                "supplied separated by ';'. The first list should contain "
                "the vertical (x) boundaries and the second list gives the "
                "horizontal (y) ones."
            )
            (
                "psf.bicubic.pixrej",
                opt::value<double>(),
                "A number defining individual pixels to exclude from "
                "piceewise bicubic PSF/PRF fits. Ignored for elliptical "
                "Gaussian PSF/PRF fits. Pixels with fitting residuals "
                "(normalized by the standard deviation) bigger than this "
                "value are excluded. If zero, no pixels are rejected."
            )
            (
                "psf.bicubic.initial-aperture",
                opt::value<double>(),
                "This aperture is used to derive an initial guess for the "
                "amplitudes of sources when fitting for a piecewise bicubic "
                "PSF/PRF model by doing aperture photometry assuming a "
                "perfectly flat PSF/PRF."
            )
            (
                "psf.bicubic.max-abs-amplitude-change",
                opt::value<double>(),
                "The absolute root of sum squares tolerance of the source "
                "amplitude changes in order to declare the piecewise bicubic"
                " PSF/PRF fitting converged."
            )
            (
                "psf.bicubic.max-rel-amplitude-change",
                opt::value<double>(),
                "The relative root of sum squares tolerance of the source "
                "amplitude changes in order to declare the piecewise bicubic"
                " PSF/PRF fitting converged."
            )
            (
                "psf.bicubic.smoothing",
                opt::value<double>()->default_value(Core::NaN),
                "How much smoothing penalty to impose when fitting the "
                "PSF/PRF. Omit for no smoothing. Value can be both positive "
                "and negative and will always result in smoothing (less for "
                "negative values)."
            );
        opt::options_description generic_options(
            "Options relevant for more than one component."
        );
        generic_options.add_options()
            (
                "gain,g",
                opt::value<double>(),
                "How many electrons are converted to 1ADU."
            )
            (
                "magnitude-1adu",
                opt::value<double>(),
                "The magnitude that corresponds to a flux of 1ADU on the "
                "input image."
            );
        opt::options_description source_options(
            "Options defining how to construct and select the sources to fit"
            " on."
        );
        source_options.add_options()
            (
                "src.cover-bicubic-grid",
                "If this option is passed and the PSF model is \"bicubic\", "
                "all pixels that at least partially overlap with the grid "
                "are assigned to the corresponding source. This option is "
                "ignored for sdk PSF models."
            )
            (
                "src.min-signal-to-noise",
                opt::value<double>(),
                "How far above the background (in units of RMS) should "
                "pixels be to still be considered part of a source. Ignored "
                "if the piecewise bibucic PSF grid is used to select source "
                "pixels (--src.cover-bicubic-grid option)."
            )
            (
                "src.max-aperture",
                opt::value<double>(),
                "If this option is specified with a positive value, pixels "
                "are assigned to sources in circular apertures (the smallest"
                " such that all pixels that pass the signal to noise cut are"
                " still assigned to the source). If an aperture larger than "
                "this value is required, an exception is thrown."
            )
            (
                "src.max-sat-frac",
                opt::value<double>(),
                "If more than this fraction of the pixels assigned to a "
                "source are saturated, the source is excluded from the fit."
            )
            (
                "src.min-pix",
                opt::value<unsigned>(),
                "The minimum number of pixels that must be assigned to a "
                "source in order to include the source is the PSF fit."
            )
            (
                "src.max-pix",
                opt::value<unsigned>(),
                "The maximum number of pixels that car be assigned to a "
                "source before excluding the source from the PSF fit."
            )
            (
                "src.max-count",
                opt::value<unsigned>(),
                "The maximum number of sources to include in the fit for the"
                " PSF shape. The rest of the sources get their amplitudes "
                "fit and are used to determine the overlaps. Sources are "
                "ranked according to the sum of "
                "(background excess)/(pixel variance+background variance) "
                "of their individual non-saturated pixels."
            );
        opt::options_description background_options(
            "Options defining how to measure the background behing the "
            "sources."
        );
        background_options.add_options()
            (
                "bg.zero",
                opt::value<bool>(),
                "If this option is passed, the input image is assumed to "
                "alredy by background corrected."
            )
            (
                "bg.annulus,b",
                opt::value<Background::Annulus>(),
                "Specifies that an annulus with the given inner radius "
                "centered around the source should be used to estimate the "
                "background and its error."
            )
            (
                "bg.min-pix",
                opt::value<unsigned>(),
                "If a source's background is based on less than this many "
                "pixels, the source is excluded from the fit."
            );
        _cmdline_config.add(io_options)
                       .add(generic_options)
                       .add(psffit_options)
                       .add(source_options)
                       .add(background_options);
    }

    std::string Config::usage_help(const std::string &prog_name) const
    {
        return "Usage: " + prog_name + " [OPTION ...]";
    }

    void Config::check_and_finalize()
    {
        if(count("help") == 0) return;
        if(
            operator[]("psf.model").as<PSF::ModelType>() == PSF::SDK
            &&
            (
                operator[]("psf.sdk.minS").as<double>() < 0
                ||
                operator[]("psf.sdk.maxS").as<double>()
                <
                operator[]("psf.sdk.minS").as<double>()
            )
        )
            throw Error::CommandLine("Invalid S range for SDK fitting "
                                     "(need max(S) > min(S) > 0).");
        if(operator[]("psf.sdk.fit-tolerance").as<double>() <= 0)
            throw Error::CommandLine("Fit tolerance must be strictly "
                                     "positive.");
        if(operator[]("psf.max-chi2").as<double>() <= 0)
            throw Error::CommandLine("Maximum source chi squared must be "
                                     "strictly positive.");
        if(operator[]("psf.bicubic.pixrej").as<double>() < 0)
            throw Error::CommandLine("Pixel rejection parameter must be "
                                     "non-negative.");
        if(operator[]("gain").as<double>() <= 0)
            throw Error::CommandLine("The gain must be strictly positive.");
        if(operator[]("src.min-signal-to-noise").as<double>() <= 0)
            throw Error::CommandLine("Minimum signal to noise must be "
                                     "strictly positive.");
        if(operator[]("src.max-sat-frac").as<double>() < 0)
            throw Error::CommandLine("Maximum saturation fraction must not "
                                     "be negative.");
        if(operator[]("src.max-count").as<unsigned>() == 0)
            throw Error::CommandLine("Source limit must be strictly "
                                     "positive.");
        if(operator[]("src.max-pix").as<unsigned>() <= 0)
            throw Error::CommandLine("Maximum pixels per source must be "
                                     "strictly positive.");
        if(operator[]("bg.min-pix").as<unsigned>() <= 1)
            throw Error::CommandLine("Background pixel limit must be "
                                     "greater than 1.");
        if(
            operator[]("psf.max-iterations").as<int>() == 0
            &&
            (*this).count("io.initial-guess") == 0
        )
            throw Error::CommandLine("Requested only amplitude fitting with "
                                     "no initial guess for the PSF map.");

        PSF::EllipticalGaussian::set_default_precision(
            operator[]("psf.sdk.rel-int-precision").as<double>(),
            operator[]("psf.sdk.abs-int-precision").as<double>()
        );
        PSF::EllipticalGaussian::set_default_max_exp_coef(
            operator[]("psf.sdk.max-exp-coef").as<double>()
        );
    }
}
