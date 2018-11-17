/**\file
 *
 * \brief Define a class used when fitting for PSF amplitude with saturated
 * pixels.
 *
 * \ingroup FitPSF
 */

#ifndef __AMPLITUDE_SATURATED_PIXEL_H
#define __AMPLITUDE_SATURATED_PIXEL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../PSF/EllipticalGaussianIntegralByOrder.h"
#include <valarray>

namespace FitPSF {

    /**\brief A structure to hold the information about saturated pixels.
     *
     * Required when fitting for the amplitude of a source
     *
     * \ingroup FitPSF
     */
    class LIB_LOCAL AmplitudeSaturatedPixel {
    private:
        double __measured; ///<The measured value of the saturated pixel

        double __variance; ///<The estimated variance of the saturated pixel

        ///\brief The integral of the PSF over the saturated pixel and its
        ///derivatives (if requested).
        std::valarray<double> __psf_integral;
    public:
        ///\brief Create a saturated pixel object.
        AmplitudeSaturatedPixel(
                ///The measured value of the pixel.
                double measured,

                ///The estimated variance of the pixel.
                double variance,

                ///the PSF integral (and its derivatives) over the pixel.
                std::valarray<double> psf_integral
        ) :
            __measured(measured),
            __variance(variance),
            __psf_integral(psf_integral)
        {}

        ///The measured value of the saturated pixel
        inline double measured() const {return __measured;}

        ///The estimated variance of the saturated pixel
        inline double variance() const {return __variance;}

        ///\brief The integral of the PSF over the saturated pixel and all of
        ///its derivatives.
        inline const std::valarray<double> &psf_integral() const
        {return __psf_integral;}

        ///\brief The integral of the PSF over the saturated pixel or one of
        ///its derivatives
        inline double psf_integral(
            ///The derivate to return.
            PSF::SDKDerivative deriv
        ) const
        {return __psf_integral[deriv];}

        ///The amplitude below which this pixel should not be saturated
        double critical_amplitude() const
        {return __measured / __psf_integral[PSF::NO_DERIV];}

        ///\brief Comparison of saturated pixels based on the amplitude below
        ///which they should not be saturated
        bool operator<(
            ///The saturated pixel to compare to.
            const AmplitudeSaturatedPixel &rhs
        ) const
        {return critical_amplitude() < rhs.critical_amplitude();}
    }; //End AmplitudeSaturatedPixel class.

}//End FitPSF namespace.

#endif
