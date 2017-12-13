/**\file
 * \brief The implementation of the FitSubpix tool.
 *
 * \ingroup FitSubpix
 */

#include "FitSubpix.h"
#include "ChiSquared.h"
#include "MultiNest_fitting.h"
#include "GSL_fitting.h"
#include "NR_fitting.h"
#include "Error.h"
#include "SubPixPhotIO.h"
#include <iostream>

void FitSubPixConfig::describe_options()
{
    _hidden.add_options()
        (
         "io.image-series",
         opt::value<StringList>(),
         "A comma separated list of images containing the same "
         "(approximately) sources, whose flux should not vary "
         "(approximately) from image to image."
        )
        ("io.psffit-series",
         opt::value<StringList>(),
         "The psf fits of the input image series files. Each file should "
         "contain a list of the sources with positions and backgrounds "
         "defined, and a PSF map."
        );
    _positional.add("io.images-series", 1);
    _positional.add("io.psffit-series", 1);
    opt::options_description io_options("Input/Output options");
    io_options.add_options()
        (
         "io.output,o",
         opt::value<std::string>(),
         "Filename for the best fit sub-pixel sensitivity map."
        );
    opt::options_description fit_subpix_options(
        "Options defining the subpixel structure fit."
    );
    fit_subpix_options.add_options()
        (
         "subpix.x-split,x",
         opt::value<unsigned>(),
         "The number of pieces in which to split pixels in the x direction."
        )
        (
         "subpix.y-split,y",
         opt::value<unsigned>(),
         "The number of pieces in which to split pixels in the y direction."
        )
    _cmdline_config.add(io_options)
                   .add(fit_subpix_options);
}

std::string FitSubPixConfig::usage_help(const std::string &prog_name) const
{
    return ("Usage: "
            + 
            prog_name
            +
            " [OPTION ...] <FITS>,<FITS>[,<FITS>...] "
            "<psf fit>,<psf fit>[,<psf fit>...]");
}

void FitSubPixConfig::check_consistency()
{
    std::ostringstream msg;
    if(count("io.image-series")!=count("io.psffit-series")) {
        msg << "Unequal number of input images (" << count("io.image-series")
            << ") and psf fit files (" << count("io.psffit-series")
            << ") specified!";
        throw Error::CommandLine(msg.str());
    }
    if(count("io.image-series")<2) {
        msg << "At least 2 input images must be used to fit for the "
            " sub-pixel structure, " << count("io.image-series")
            << " found!";
        throw Error::CommandLine(msg.str());
    }
    if(
            (operator[]("subpix.x-split").as<unsigned>()*
             operator[]("subpix.y-split").as<unsigned>())
            <
            2
    ) {
        msg << "The sub-pixel map must contain at least two pieces, "
            << operator[]("subpix.x-split").as<unsigned>()
            << "x"
            << operator[]("subpix.y-split").as<unsigned>()
            << " map specified!";
        throw Error::CommandLine(msg.str());
    }
}

///Fits for the subpixel structure in a set of images.
int main(int argc, char **argv)
{
	try {
		FitSubPixConfig options(argc, argv);
		if(!options.proceed()) return 1;

        std::string fit_method=options["subpix.method"].as<std::string>();
		if(fit_method=="simplex")
			fit_using_GSL_simplex(options);
		else if(fit_method=="siman")
			fit_using_GSL_simulated_annealing(options);
		else if(fit_method=="multinest")
			MultiNestFit::fit(options["subpix.x-split"].as<unsigned>(), 
                              options["subpix.y_split"].as<unsigned>(),
                              options["io.image-series"].as<StringList>(),
                              options["io.psffit-series"].as<StringList>(),
                              options["subpix.aperture"].as<double>());
		else fit_using_NR(options);
	} catch(Error::General &ex) {
		std::cerr << ex.what() << ":" << ex.get_message() << std::endl;
	}
}
