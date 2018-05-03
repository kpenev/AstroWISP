/**\file
 *
 * \brief Declares a base class for sources participating in SDK fitting.
 *
 * \ingroup FitPSF
 */

#ifndef __SDK_SOURCE_BASE_H
#define __SDK_SOURCE_BASE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "AmplitudeSaturatedPixel.h"
#include "Source.h"
#include "../PSF/EllipticalGaussianIntegralByOrder.h"
#include "../Background/Source.h"
#include "../Core/Image.h"
#include "../Core/SourceID.h"
#include <list>
#include <string>

namespace FitPSF {

    /**\brief An abstract class that defines the minimum requirements of a 
     * source to participate in the fitting for an elliptical gaussian PSF 
     * model.
     *
     * \ingroup FitPSF
     */
    class LIB_LOCAL SDKSourceBase : public Source< PSF::EllipticalGaussian > {
    private:
        ///\brief Should first order S,D,K derivatives of the amplitude and
        /// \f$\chi^2\f$ be calculated
        bool __first_deriv, 

             ///\brief Should second order S,D,K derivatives of the amplitude 
             ///and \f$\chi^2\f$ be calculated
             __second_deriv;
             
        std::valarray<double>
            ///The best fit amplitude (and derivatives) for the source PSF
            __amplitude,
            
            ///\brief The \f$\chi^2\f$ residual of the best amplitude fit and 
            ///its derivatives
            __chi2;

        ///\brief Adds a pixel to the source.
        ///
        ///Adds to the numerator and denominator appropriate quantities
        ///calculated from either the current pixel (if sat_iter is 
        ///__saturated.rend()), or the given saturated entry if it is.
        void add_pixel_to_scaling_fit(
            const PSF::EllipticalGaussian &psf,
            std::valarray<double> &numerator,
            std::valarray<double> &denominator,
            double &measured2_term,
            const std::list<AmplitudeSaturatedPixel>::const_reverse_iterator
            &saturated_iter
        );

    protected:
        ///\brief Calculates the values and, if requested, derivatives of
        /// \f$\chi^2\f$ and the amplitude from the given quantities.
        virtual void fit_scaling(
            double measured2_term,
            const std::valarray<double> &numerator,
            const std::valarray<double> &denominator
        );

    public:
        ///Create a source for fitting an elliptical gaussian PSF model.
        template<class EIGEN_MATRIX>
            SDKSourceBase(
                ///See Source::Source()
                const Core::SourceID &id,

                ///See Source::Source()
                double x0,

                ///See Source::Source()
                double y0,

                ///The background under the pixel (in ADU)
                const Background::Source &background,

                ///The actual image we are deriving a PSF map for
                const Core::Image<double> &observed_image,

                ///A two dimensional array which keeps track of what pixels 
                ///of the input image are assigned to what source. On exit it 
                ///is updated with the pixels belonging to the newly 
                ///constructed source.
                EIGEN_MATRIX &source_assignment,

                ///The gain (electrons per ADU) in the observed image
                double gain,

                ///How much above the background a pixel needs to be in order 
                ///to be allocated to this source (the alpha parameter in the
                ///description)
                double alpha,

                ///The id to assign to this source in the source_assignment 
                ///array
                int source_id,

                ///See Source constructor.
                const std::string &output_fname,

                ///Should first order derivatives with respect to S, D and K 
                ///be calculated
                bool calculate_first_deriv=false,

                ///Should second order derivatives with respect to S, D and K 
                ///be calculated
                bool calculate_second_deriv=false,

                ///If nonzero impose a circular aperture for the source no 
                ///larger than the given value (otherwise uses only pixels 
                ///inconsistent with the background at the prescribed by 
                ///alpha level). The size of the circular aperture is the 
                ///smallest size possible that encapsulates all pixels that 
                ///pass the alpha test.
                double max_circular_aperture=0
            ) : 
                Source(id, x0, y0,
                       background,
                       observed_image,
                       source_assignment,
                       gain,
                       alpha,
                       source_id,
                       max_circular_aperture,
                       output_fname),
                __first_deriv(calculate_first_deriv),
                __second_deriv(calculate_second_deriv),
                __amplitude(PSF::KK_DERIV+1),
                __chi2(PSF::KK_DERIV+1)
        {}

        ///\brief Enables the calculation of the first order derivatives of 
        ///the amplitude and \f$\chi^2\f$.
        void enable_first_deriv() {__first_deriv = true;}

        ///\brief Disables the calculation of the first order derivatives of 
        ///the amplitude and \f$\chi^2\f$.
        void disable_first_deriv() {__first_deriv = false;}

        ///\brief Enables the calculation of the second order derivatives of 
        ///the amplitude and \f$\chi^2\f$.
        void enable_second_deriv() {__second_deriv = true;}

        ///\brief Disables the calculation of the second order derivatives of 
        ///the amplitude and \f$\chi^2\f$.
        void disable_second_deriv() {__second_deriv = false;}

        ///\brief Should first order derivatives of the amplited and 
        /// \f$\chi^2\f$ be calculated.
        bool calculate_first_deriv() const {return __first_deriv;}

        ///\brief Should second order derivatives of the amplited and
        /// \f$\chi^2\f$ be calculated.
        bool calculate_second_deriv() const {return __second_deriv;}

        ///Returns the amplitude of the last scaling fit
        double amplitude(PSF::SDKDerivative deriv = PSF::NO_DERIV)
        {return __amplitude[deriv];}

        ///Returns the \f$\chi^2\f$ of the last scaling fit or any derivative.
        double chi2(PSF::SDKDerivative deriv) const {return __chi2[deriv];}

        ///Returns the \f$\chi^2\f$ of the last scaling fit
        double chi2() const {return __chi2[PSF::NO_DERIV];}

        ///\brief Returns the \f$\chi^2\f$ of the last scaling fit and its up 
        ///to second order derivatives.
        const std::valarray<double> &chi2_all_deriv() const
        {return __chi2;}

        ///\brief Fits for the best scaling between the PSF integrals over 
        ///pixels and the measured pixel values.
        ///
        ///Initializes the __amplitude member with the best fit scaling and 
        ///the __chi2 member with the reduced \f$\chi^2\f$ of the best fit 
        ///scaling.
        void fit_scaling(const PSF::EllipticalGaussian &psf);
    };

} //End FitPSF namespace.

#endif
