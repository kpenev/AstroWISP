#ifndef __FIT_SUBPIX_COMMAND_LINE_OPTIONS_H
#define __FIT_SUBPIX_COMMAND_LINE_OPTIONS_H

/**\file
 * \brief Declares the FitSubpixCommandLineOptions parsers for the FitSubpix
 * tool.
 *
 * \ingroup FitSubpix
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "PhotColumns.h"
#include "CommandLineUtil.h"
#include "SubPixelMap.h"
#include <vector>
#include <string>
#include <list>
#include <valarray>
#include <gsl/gsl_vector.h>

///Tags for the fitting methods available.
enum LIB_LOCAL FIT_METHOD {
    MultiNest, 
    GSL_simplex, 
    NewtonRaphson,
	GSL_simulated_annealing
};

///Tags for the various simulated annealing options.
enum LIB_LOCAL SIMAN_OPTION {
    NTRIES,
    MAX_STEP,
    BOLTZMAN_K,
    START_TEMPERATURE,
	COOLING_RATE,
    MIN_TEMPERATURE
};

///\brief Everything defined through the command line is accessible as a
///method.
///
///\ingroup FitSubpix
class LIB_LOCAL FitSubpixCommandLineOptions {
private:
	///The input fits frame filenames
	std::vector<std::string> __frames;

	///The files containing lists of sources for each frame.
	std::vector<std::string> __sources;

	
	unsigned __x_split, ///< The number of x splits of a pixel
			 __y_split, ///< The number of y splits of a pixel

			 ///The number of points to try for each simulated annealing step
			 __siman_ntries;

	///The aperture to use when doing photometry.
	double __aperture,

		   ///The maximum step size in the simulated annealing random walk
		   __siman_max_step,
		   
		   ///Boltzman constant for the simulated annealing
		   __siman_boltzman_k,
		   
		   ///Starting temperature for the simulated annealing
		   __siman_start_temperature,
		   
		   ///Simulated annealing cooling factor
		   __siman_cooling_rate,
		   
		   ///Simulated annealing minimum temperature
		   __siman_min_temperature,
		   
		   ///The minimum sudbivision to impose when integrating subpixels
		   __max_exp_coef;

	///The annulus around a source on which to base the background estimate
	BackgroundAnnulus __bg_annulus;

	///The order of the columns in the input
	std::list<Phot::Columns> __input_columns;

	///What method to use for fitting for the subpixel structure.
	FIT_METHOD __fit_method;
	
	///\brief True only if no problems were encoutered while parsing the
	///command line.
	bool __parsed_ok;
public:
	///Parse the command line.
	FitSubpixCommandLineOptions(int argc, char **argv);

	///The fits frames to do photometry on.
	const std::vector<std::string> &frames() const {return __frames;}

	///The sources on the fits frames to do photometry for.
	const std::vector<std::string> &sources() const {return __sources;}

	///In how many parts to split the pixel in the x direction
	unsigned x_split() const {return __x_split;}

	///In how many parts to split the pixel in the y direction
	unsigned y_split() const {return __y_split;}

	///What aperture to use when doing photometry.
	double aperture() const {return __aperture;}

	///The minimum sudbivision to impose when integrating subpixels
	double max_exp_coef() const {return __max_exp_coef;}

	///The annulus around a source on which to base the background estimate
	const BackgroundAnnulus &background_annulus() const
	{return __bg_annulus;}

	///The order of the columns in the input
	const std::list<Phot::Columns> &input_columns() const
	{return __input_columns;}

	///Returns the method that should be used for fitting for the subpixel 
	///structure.
	FIT_METHOD fit_method() {return __fit_method;}

	///The value of a single simulated annealing parameter to use for fitting.
	double siman_option(SIMAN_OPTION option) const;

	///Did command line parsing occur without any problems?
	operator bool() const {return __parsed_ok;}
};

#endif
