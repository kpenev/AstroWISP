/**\file
 *
 * \brief Defines some of the methods of the SDKSourceBase class.
 *
 * \ingroup FitPSF
 */

#include "SDKSourceBase.h"

namespace FitPSF {

    ///\brief Takes a second derivative and returns the two first order 
    ///derivatives which it is composed of (e.g. PSF::SS_DERIV ->
    ///PSF::S_DERIV, PSF::S_DERIV). If deriv is not a second derivative,
    ///the behaviour is undefined.
    void split_second_deriv(PSF::SDKDerivative second,
                            PSF::SDKDerivative &deriv1,
                            PSF::SDKDerivative &deriv2)
    {
        if(second == PSF::SS_DERIV) deriv1 = deriv2 = PSF::S_DERIV;
        else if(second == PSF::DD_DERIV) deriv1 = deriv2 = PSF::D_DERIV;
        else if(second == PSF::KK_DERIV) deriv1 = deriv2 = PSF::K_DERIV;
        else if(second == PSF::SD_DERIV) {
            deriv1 = PSF::S_DERIV;
            deriv2 = PSF::D_DERIV;
        } else if(second == PSF::SK_DERIV) {
            deriv1 = PSF::S_DERIV;
            deriv2 = PSF::K_DERIV;
        } else if(second == PSF::DK_DERIV) {
            deriv1 = PSF::D_DERIV;
            deriv2 = PSF::K_DERIV;
        }
        else 
            throw Error::InvalidArgument("PSFFitting.cpp:split_second_deriv",
                                         "Unrecognized second derivative.");
    }

    void SDKSourceBase::add_pixel_to_scaling_fit(
        const PSF::EllipticalGaussian &psf,
        std::valarray<double> &numerator,
        std::valarray<double> &denominator,
        double &measured2_term,
        const std::list<AmplitudeSaturatedPixel>::const_reverse_iterator
            &saturated_iter
    )
    {
        double variance = (current_pixel()->variance()
                           +
                           background_electrons_variance()),
               measured = (current_pixel()->measured()
                           -
                           background_electrons()),
               psf_no_deriv;
        measured2_term += std::pow(measured, 2) / variance;

        for(int d = PSF::NO_DERIV; d <= PSF::KK_DERIV; ++d)
            if((d <= PSF::K_DERIV && __first_deriv) || __second_deriv) {
                PSF::SDKDerivative
                    deriv = static_cast<PSF::SDKDerivative>(d);
                double psf_deriv = (
                    saturated_iter == _saturated_pixels.rend()
                    ? pixel_psf(psf, deriv)
                    : saturated_iter->psf_integral(deriv)
                );
                if(deriv == PSF::NO_DERIV) psf_no_deriv = psf_deriv;
#ifdef DEBUG
                assert(!std::isnan(psf_deriv));
#endif
                numerator[deriv] += measured * psf_deriv / variance;
                if(deriv == PSF::NO_DERIV)
                    denominator[deriv] += (std::pow(psf_no_deriv, 2)
                                           /
                                           variance);
                else if(deriv <= PSF::K_DERIV)
                    denominator[deriv] += (
                        2.0 * psf_no_deriv * psf_deriv / variance
                    );
                else {
                    PSF::SDKDerivative deriv1, deriv2;
                    split_second_deriv(deriv, deriv1, deriv2);
                    double psf_deriv1_deriv2;
                    if(saturated_iter == _saturated_pixels.rend())
                        psf_deriv1_deriv2 = (pixel_psf(psf, deriv1)
                                             *
                                             pixel_psf(psf, deriv2));
                    else psf_deriv1_deriv2 = (
                        saturated_iter->psf_integral(deriv1)
                        *
                        saturated_iter->psf_integral(deriv2)
                    );
                    denominator[deriv] += 2.0 * (
                        psf_no_deriv * psf_deriv + psf_deriv1_deriv2
                    ) / variance;
                }
            }
    }

    void SDKSourceBase::fit_scaling(
        double measured2_term,
        const std::valarray<double> &numerator,
        const std::valarray<double> &denominator
    )
    {

        flux()[0].value() = (numerator[PSF::NO_DERIV]
                             /
                             denominator[PSF::NO_DERIV]);
        __chi2[PSF::NO_DERIV] = (measured2_term
                                 -
                                 std::pow(numerator[PSF::NO_DERIV], 2)
                                 /
                                 denominator[PSF::NO_DERIV]);
        flux()[0].error() = std::sqrt(chi2()
                                   /
                                   denominator[PSF::NO_DERIV]
                                   /
                                   (flux_fit_pixel_count() - 1));
        __amplitude[PSF::NO_DERIV] = flux()[0].value();
        if(__first_deriv) 
            for(int d = PSF::S_DERIV; d <= PSF::K_DERIV; ++d) {
                PSF::SDKDerivative
                    deriv = static_cast<PSF::SDKDerivative>(d);
                __amplitude[deriv] = (
                    numerator[deriv] / denominator[PSF::NO_DERIV]
                    -
                    numerator[PSF::NO_DERIV]
                    /
                    std::pow(denominator[PSF::NO_DERIV], 2)
                    *
                    denominator[deriv]
                );
                __chi2[deriv] = (
                    -__amplitude[deriv] * numerator[PSF::NO_DERIV]
                    -
                    __amplitude[PSF::NO_DERIV] * numerator[deriv]
                );
            }
        if(__second_deriv) {
            double denom2 = std::pow(denominator[PSF::NO_DERIV], 2),
                   denom3 = denom2 * denominator[PSF::NO_DERIV];
            for(int d = PSF::SS_DERIV; d <= PSF::KK_DERIV; ++d) {
                PSF::SDKDerivative
                    deriv = static_cast<PSF::SDKDerivative>(d), 
                    deriv1,
                    deriv2;
                split_second_deriv(deriv, deriv1, deriv2);
                __amplitude[deriv] = (
                    numerator[deriv] / denominator[PSF::NO_DERIV]
                    -
                    (
                        numerator[deriv1] * denominator[deriv2]
                        +
                        numerator[deriv2] * denominator[deriv1]
                    ) / denom2
                    +
                    (
                        2.0
                        *
                        numerator[PSF::NO_DERIV]
                        *
                        denominator[deriv1]
                        *
                        denominator[deriv2]
                        /
                        denom3
                    )
                    -
                    numerator[PSF::NO_DERIV] * denominator[deriv] / denom2
                );
                __chi2[deriv] = (
                    -__amplitude[deriv1] * numerator[deriv2]
                    -
                    __amplitude[deriv2] * numerator[deriv1]
                    -
                    __amplitude[deriv] * numerator[PSF::NO_DERIV]
                    -
                    __amplitude[PSF::NO_DERIV] * numerator[deriv]
                );
#ifdef DEBUG
                assert(!std::isnan(__amplitude[deriv]));
                assert(!std::isnan(__chi2[deriv]));
#endif
            }
        }
    }

    ///\todo Re-enable saturated pixel handling.
    void SDKSourceBase::fit_scaling(const PSF::EllipticalGaussian &psf)
    {
        double measured2_term = 0.0; 
        std::valarray<double> numerator(0.0, PSF::KK_DERIV + 1),
                              denominator(0.0, PSF::KK_DERIV + 1);
        fit_pixels() = 0;
        _saturated_pixels.clear();
        restart_pixel_iteration();
        do {
            if(current_pixel_flag() == Core::GOOD) {
                add_pixel_to_scaling_fit(psf,
                                         numerator,
                                         denominator,
                                         measured2_term,
                                         _saturated_pixels.rend());
                ++fit_pixels();
            } else {
    /*			std::valarray<double> psf_integral(Core::NaN, PSF::KK_DERIV+1);
                psf_integral[PSF::NO_DERIV]=pixel_psf();
                for(int d=(__first_deriv ? PSF::S_DERIV : PSF::SS_DERIV);
                        d<=(__second_deriv ? PSF::KK_DERIV : PSF::K_DERIV); d++)
                    psf_integral[d]=pixel_psf(static_cast<PSF::SDKDerivative>(d));
                _saturated_pixels.push_back(
                        AmplitudeSaturatedPixel(pixel_measured()-
                        background_electrons(),
                            psf_integral));*/
            }
        } while(next_pixel());
        _saturated_pixels.sort();
        __amplitude[PSF::NO_DERIV] = (numerator[PSF::NO_DERIV]
                                      /
                                      denominator[PSF::NO_DERIV]);

        for(
            std::list<AmplitudeSaturatedPixel>::const_reverse_iterator
                saturated_iter = _saturated_pixels.rbegin();
            (
                saturated_iter != _saturated_pixels.rend()
                &&
                (
                    saturated_iter->critical_amplitude() 
                    >
                    __amplitude[PSF::NO_DERIV]
                )
            );
            ++saturated_iter
        ) {
            add_pixel_to_scaling_fit(psf,
                                     numerator,
                                     denominator,
                                     measured2_term,
                                     saturated_iter);
            __amplitude[PSF::NO_DERIV] = (numerator[PSF::NO_DERIV]
                                          /
                                          denominator[PSF::NO_DERIV]);
            ++fit_pixels();
        }
        fit_scaling(measured2_term, numerator, denominator);
    }

} //End FitPSF namespace.
