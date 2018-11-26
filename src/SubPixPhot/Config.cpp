/**\file
 *
 * \brief Implement the methods of the Config class.
 *
 * \ingroup SubPixPhot
 */

#include "Config.h"

namespace SubPixPhot {
    void Config::describe_options()
    {
        _hidden.add_options()
            (
                "io.image",
                opt::value<std::string>(),
                "Reduced file to fit the PSF of."
            )
            (
                "io.sources,l",
                opt::value<std::string>(),
                "HDF5 file listing the sources to perform photometry for. "
                "Typically produced by FitPSF."
            );
        _positional.add("io.image", 1);
        _positional.add("io.sources", 1);
        opt::options_description io_options("Input/Output options");
        io_options.add_options()
            (
                "io.psfmap",
                opt::value<std::string>(),
                "HDF5 file defining the PSF map. If not specified, defaults "
                "to the input sources file."
            )
            (
                "io.output,o",
                opt::value<std::string>(),
                "File to place the output in HDF5 format. If not specified, "
                "defaults to the input sources file."
            )
            (
                "io.subpix,s",
                opt::value<std::string>()->default_value(""),
                "FITS file representing the sensitivity structure of a "
                "pixel. If not specified, uniform sensitivity is assumed."
            )
            (
                "io.output-quantities",
                opt::value<Core::ColumnList>(),
                "A comma separated list of quantities to output. Recognized "
                "values: id, x, y, S, D, K, A|amp, flux, flux_err, mag, "
                "mag_err, flag, bg, bg_err, enabled, sn, npix, nbgpix"
#ifdef DEBUG
                ", time"
#endif
                ". In text output a group of consecutive columns that depend"
                " on aperture are repeated automatically."
            );
        opt::options_description apphot_options(
            "Options configuring how aperture photometry is performed.");
        apphot_options.add_options()
            (
                "ap.aperture",
                opt::value<Core::RealList>(),
                "Comma separated list of apertures to use."
            )
            (
                "ap.const-error",
                opt::value<double>(),
                "A value to add to the error estimate of a pixel (intended "
                "to represent things like readout noise, truncation noise "
                "etc.). "
            );
        opt::options_description generic_options(
            "Options relevant for more than one component.");
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
        opt::options_description background_options(
            "Options defining how to measure the background behing the "
            "sources.");
        background_options.add_options()
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
            .add(apphot_options)
            .add(background_options);
    }

    void SubPixPhot::Config::apply_fallbacks()
    {
        std::stringstream fallbacks;
        fallbacks
            << "[io]" << std::endl
            << "psfmap = " << (*this)["io.sources"].as<std::string>()
            << std::endl
            << "output = " << (*this)["io.sources"].as<std::string>()
            << std::endl;
        fallbacks.seekg(0);
        opt::store(
            opt::parse_config_file(fallbacks,
                                   _cmdline_config),
            *this
        );
        opt::notify(*this);
    }

    void SubPixPhot::Config::check_consistency()
    {
        if((*this)["ap.aperture"].as<Core::RealList>().size() == 0)
            throw Error::CommandLine("No apertures defined!");
        for(
                Core::RealList::const_iterator
                    i = (*this)["ap.aperture"].as<Core::RealList>().begin();
                i != (*this)["ap.aperture"].as<Core::RealList>().end();
                ++i
        )
            if(*i < 0)
                throw Error::CommandLine("Negative aperture encountered!");
        if((*this)["gain"].as<double>()<0)
            throw Error::CommandLine("Negative gain specified!");
        assert(count("io.psfmap"));
        assert(count("io.output"));
    }

    std::string SubPixPhot::Config::usage_help(
        const std::string &prog_name
    ) const
    {
        return "Usage: "
               +
               prog_name
               +
               " [OPTION ...] <FITS image> <HDF5 sources>";
    }

} //End SubPixPhot namespace.
