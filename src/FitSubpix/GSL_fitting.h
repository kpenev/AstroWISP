/**\file
 *
 * \brief Declarations for simplex and simulated annealing fitting for the
 * sub-pixel structure.
 *
 * \ingroup FitSubpix
 */

#ifndef __GSL_FITTING_H
#define __GSL_FITTING_H

#include "SubPixelMap.h"
#include "FittingUtil.h"
#include "FitSubpix.h"
#include "Error.h"
#include "ChiSquared.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_siman.h>
#include <iostream>
#include <iomanip>

///A subpixel map that sets its sensitivities from a GSL vector.
class GSLSubPixelMap : public SubPixelMap {
public:
	///Construct a map with the given resolution.
	GSLSubPixelMap(unsigned long x_split, unsigned long y_split) :
		SubPixelMap(x_split, y_split, "GSL sub pixel map") {}

	///\brief Sets the sensitivities from the given GSL vector.
	///
	///The convention is that the vector contains directly all sensitivity
	///values except the maximum x and y, which is set to (number of
	///subpixels - sum of all other sensitivities).
	void set_sensitivities(const gsl_vector *cube);
};

///\brief Fit for the subpixel structure using the simplex method.
///
///Minimizes the variance of a list of  sources over a set of frames by
///varying the subpixel structure used when extracting the fluxes.
void fit_using_GSL_simplex(const FitSubPixConfig &options);

///\brief Fit for the supixel structure using simulated annealing.
///
///Minimizes the variance of a list of sources over a set of frames by
///varying the subpixel structure used when extracting the fluxes.
void fit_using_GSL_simulated_annealing(const FitSubPixConfig &options);

///Outputs the subpixel map currently in the given minimizer.
std::ostream &operator<<(std::ostream &os, 
				gsl_multimin_fminimizer *minimizer);

#endif
