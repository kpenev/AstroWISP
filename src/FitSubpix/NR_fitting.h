/**\file
 *
 * \brief Declares the function that performs Newton-Raphson fitting for the
 * sub-pixel sensitivity map.
 *
 * \ingroup FitSubpix
 */

#ifndef __NR_FITTING_H
#define __NR_FITTING_H

#include "../Core/SharedLibraryExportMacros.h"
#include "FittingUtil.h"
#include "FitSubpixCommandLineOptions.h"
#include "NaN.h"

///\brief Uses Newton-Raphson method to find a zero in the derivative of the
///variance w.r.t. the sub-pixel sensitivity map.
LIB_LOCAL void fit_using_NR(const FitSubpixCommandLineOptions &options);

#endif
