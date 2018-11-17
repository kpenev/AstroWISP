/**\file
 *
 * \brief Define the functions from SDKUtil.h.
 *
 * \ingroup FitPSF
 */

#include "SDKUtil.h"
#include <sstream>

namespace FitPSF {

    template<class SOURCE_ITERATOR>
        double calculate_chi2(SOURCE_ITERATOR first_source,
                              SOURCE_ITERATOR past_last_source,
                              double max_source_chi2)
        {
            double numerator = 0,
                   denominator = 0;
            for(
                SOURCE_ITERATOR si = first_source;
                si!=past_last_source; ++si
            ) {
                if(std::isnan(si->chi2())) {
                    std::ostringstream msg;
                    msg << "NaN chi2 for " << si->id()
                        << ", bg=" << si->background_electrons()
                        << ", bg variance="
                        << si->background_electrons_variance()
                        << ", "
                        << si->pixel_count() << " assigned pixels.";
                    throw Error::Fitting(msg.str());
                }
                if(si->chi2() / (si->fit_pixels() - 1) <= max_source_chi2) {
                    numerator += si->chi2();
                    denominator += si->fit_pixels() - 1;
                }
            }
            return numerator / denominator;
        }

    double gsl_minimization_function(const gsl_vector *poly_coef,
                                     void *params)
    {
        void **param_array = reinterpret_cast<void**>(params);
        GSLSourceIteratorType 
            first_source = *reinterpret_cast<GSLSourceIteratorType*>(
                    param_array[0]
            ),
            past_last_source = *reinterpret_cast<GSLSourceIteratorType*>(
                param_array[1]
            );
        PolynomialSDK<GSLSubPixType> 
            *poly_sdk = reinterpret_cast< PolynomialSDK<GSLSubPixType>* >(
                param_array[2]
            );
        const GSLSubPixType *subpix_map = reinterpret_cast<GSLSubPixType*>(
            param_array[3]
        );
        double max_source_chi2 = *reinterpret_cast<double*>(param_array[4]),
               minS = *reinterpret_cast<double*>(param_array[5]),
               maxS = *reinterpret_cast<double*>(param_array[6]);
        std::valarray<double> poly_coef_array(poly_coef->size);
        for(size_t i = 0; i < poly_coef->size; ++i)
            poly_coef_array[i] = gsl_vector_get(poly_coef, i);
        for(
            GSLSourceIteratorType src_iter = first_source;
            src_iter != past_last_source;
            ++src_iter
        ) {
            PSF::EllipticalGaussian psf = (*poly_sdk)(*src_iter,
                                                      poly_coef_array);
            if(psf.s() < minS || psf.s() > maxS)
                return std::numeric_limits<double>::max();
            src_iter->set_PSF(psf, *subpix_map);
        }
        return calculate_chi2(first_source,
                              past_last_source,
                              max_source_chi2);
    }

    double gsl_s_minimization_function(double s, void *params)
    {
        void **param_array = reinterpret_cast<void**>(params);
        GSLSourceIteratorType 
            first_source = *reinterpret_cast<GSLSourceIteratorType*>(
                param_array[0]
            ),
            past_last_source = *reinterpret_cast<GSLSourceIteratorType*>(
                param_array[1]
            );
        const GSLSubPixType *subpix_map = reinterpret_cast<GSLSubPixType*>(
            param_array[2]
        );
        PSF::EllipticalGaussian psf(s, 0.0, 0.0, 0, 1e-5, 1e-3);
        for(
            GSLSourceIteratorType src_iter = first_source;
            src_iter != past_last_source;
            ++src_iter
        ) {
            src_iter->set_PSF(psf, *subpix_map);
        }
        return calculate_chi2(first_source, past_last_source);
    }

} //End FitPSF namespace.
