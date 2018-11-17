/**\file
 *
 * \brief Declare functions used during elliptical Gaussian PSF fitting.
 *
 * \ingroup FitPSF
 */

#ifndef __SDK_UTIL_H
#define __SDK_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SDKSource.h"
#include <gsl/gsl_vector.h>

namespace FitPSF {

    template< class SUBPIX_TYPE > class PolynomialSDK;

    ///The \f$\chi^2\f$ that corresponds to a sequence of sources.
    template<class SOURCE_ITERATOR>
        double LIB_LOCAL calculate_chi2(
            ///A source iterator pointing to the first source in the
            ///sequence.
            SOURCE_ITERATOR first_source,

            ///A source iterator pointing to immediately past the last source 
            ///in the sequence.
            SOURCE_ITERATOR past_last_source,

            ///If \f$\chi^2\f$ exceeds this value, the source is considered
            ///non-point.
            double max_source_chi2=Core::Inf
        );

    ///\brief The function to pass to the GSL simplex minimizer, forced to 
    ///assume SOURCE_ITERATOR = GSLSourceIteratorType
    ///and SUBPIX_TYPE = GSLSubPixType
    ///
    ///See GSL documentation for description of the parameters.
    LIB_LOCAL double gsl_minimization_function(
        const gsl_vector *poly_coef,
        void *params
    );

    ///\brief The function to pass to the GSL brent minimizer for finding an 
    ///initial frame S.
    ///
    ///See GSL documentation for description of the parameters.
    LIB_LOCAL double gsl_s_minimization_function(double s, void *params);

} //End FitPSF namespace.

#endif
