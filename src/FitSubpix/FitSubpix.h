#ifndef __FIT_SUBPIX_H
#define __FIT_SUBPIX_H

/**\file
 * \brief Declares the FitSubpixCommandLineOptions parsers for the FitSubpix
 * tool.
 *
 * \ingroup FitSubpix
 */

#include "PhotColumns.h"
#include "CommandLineUtil.h"
#include "SubPixelMap.h"
#include <vector>
#include <string>
#include <list>
#include <valarray>
#include <gsl/gsl_vector.h>

///Default configuration from file but overwritten by command line options.
class FitSubPixConfig : public CommandLineConfig {
private:
	///Describes the available command line options.
	void describe_options();

	///The part of the help describing the usage and purpose (no options).
	std::string usage_help(const std::string &prog_name) const;

	///\brief Checks for consistency between the command line options.
	///
	///Throws an exception if some inconsistency is detected.
	void check_consistency();
public:
	///Parse the command line.
	FitSubPixConfig(
			///The number of arguments on the command line
			///(+1 for the executable)
			int argc,
			
			///A C style array of the actual command line arguments.
			char **argv)
	{
		parse(argc, argv);
		if(count("help")==0) {
			check_consistency();
		}
        PSF::EllipticalGaussian::set_default_precision(
            operator[]("psf.sdk.rel-int-precision").as<double>(),
            operator[]("psf.sdk.abs-int-precision").as<double>()
        );
        PSF::EllipticalGaussian::set_default_max_exp_coef(
            operator[]("psf.sdk.max-exp-coef").as<double>()
        );
	}
};

#endif
