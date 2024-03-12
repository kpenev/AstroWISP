/**\file
 *
 * \brief Defines the base class for sources participating in PSF fitting.
 *
 * \ingroup FitPSF
 */

#ifndef __FITPSF_SOURCE_H
#define __FITPSF_SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "AmplitudeSaturatedPixel.h"
#include "Common.h"
#include "../PSF/EllipticalGaussian.h"
#include "../PSF/Map.h"
#include "../PSF/MapSource.h"
#include "../Background/Source.h"
#include "../Core/Source.h"
#include "../Core/SourceID.h"
#include "../Core/Flux.h"
#include "../Core/Point.h"
#include "../Core/SubPixelMap.h"
#include "../Core/Typedefs.h"
#include "../IO/H5IODataTree.h"
//#include "SubPixPhotIO.h"
#include "Eigen/Dense"
#include <list>
#include <set>
#include <valarray>
#include <iostream>
#include <iomanip>

namespace FitPSF {

    /**\brief An abstract base class for sources which participate in PSF
     * fitting.
     *
     * \ingroup FitPSF
     */
    template<class PSF_TYPE>
        class LIB_LOCAL Source : public PSF::MapSource {
        private:

            double
                ///The background under the source in electrons
                __background_electrons,

                ///The variance of the background under the source
                __background_electrons_variance,

                ///The gain to assume for image pixels
                __gain,

               ///The best fit amplitude of this source
               __amplitude;


            ///The id of the source in the source_assignment array
            unsigned __source_assignment_id,

                     ///\brief The number of sources which were combined with this
                     ///source when fitting for the amplitude.
                     __group_sources,

                     ///Identifier of the image this source is part of
                     __image_id;

            ///See subpix_map argument of the constructors.
            const Core::SubPixelMap *__subpix_map;

            ///Is this source rejcted as non-point source.
            bool __nonpoint;

            ///\brief Source flux estimated by scaling for flux not in the source
            ///pixels.
            Core::Flux __mask_flux;

            std::string
                ///See output_filename().
                __output_filename;

            ///\brief The terms participating in the expansion of each of the PSF
            ///parameters.
            Eigen::VectorXd __expansion_terms;

            ///See chi2() method.
            double __chi2;

            ///The reason the source was dropped if it was.
            SourceDropReason __drop_reason;

        protected:
            ///\brief A structure holding all necessary information about
            ///saturated pixels.
            std::list<AmplitudeSaturatedPixel> _saturated_pixels;

            ///\brief Calculate the dot of pixel excesses with the given quantity
            ///as well as with itself.
            ///
            ///Having this function virtual allows optimizing the calculation in
            ///some cases (e.g. when fitting for PSF which linearly depends on
            ///its parameters).
            virtual void pixel_excess_reductions(
                ///A vector of values to dot with the pixel excess values.
                const Eigen::VectorXd &to_dot_with,

                ///On exit this variable is set to the dot product of to_dot_with
                ///and the pixel excesses.
                double &dot_product,

                ///On exit this variable is set to the sum of the squares of the
                ///pixel excesses.
                double &excess_sum_squares
            ) =0;

            ///\brief Fit for the flux of the source and return the change from
            ///its previous value.
            double fit_flux(
                ///The integrals of the PSF x sub-pixel map over the source
                ///pixels used for flux fitting divided by the variance.
                const Eigen::VectorXd &estimated_excesses
            );

        public:
            ///Creates a PSF fitting source.
            Source(
                ///The ID of the source.
                const Core::SourceID &id,

                ///The x coordinate of the source center.
                double x0,

                ///The y coordinate of the source center.
                double y0,

                ///The gain to assume for the backcground measurement.
                double gain,

                ///The background under the source (in ADU)
                const Background::Source &background,

                ///The id to assign to this source in the source_assignment array
                int source_assignment_id,

                ///The sub-pixel sensitivity map to assume. Must not be destroyed
                ///while this object is in use.
                const Core::SubPixelMap *subpix_map,

                ///The file name of the image from which this source was
                ///extracted.
                unsigned image_id,

                ///The name of the file where this source should be saved after
                ///the fit.
                const std::string &output_fname
            );

            ///Is this source ready for fitting?
            virtual bool ready_to_fit() const {return true;}

            ///The background under the source in electrons
            double background_electrons() const {return __background_electrons;}

            ///Set the background under the source in electrons
            void set_background_electrons(double value)
            {
                __background_electrons=value;
                background().value() = value / __gain;
            }

            ///\brief The variance of the background under the source in
            ///electrons squared
            double background_electrons_variance() const
            {return __background_electrons_variance;}

            ///\brief A reference to the variance of the background under the
            ///source in electrons squared
            void set_background_electrons_variance(double value)
            {
                __background_electrons_variance = value;
                background().error() = std::sqrt(value) / __gain;
            }

            ///The number of pixels used when measuring the source background.
            unsigned background_pixels() const {return background().pixels();}

            ///\brief A reference to the number of pixels used when measuring the
            ///source background.
            unsigned &background_pixels() {return background().pixels();}

            ///Calculate and return the __mask_flux using the given PSF/PRF.
            virtual const Core::Flux &calculate_mask_flux(const PSF::PSF &psf) = 0;

            ///See __mask_flux (filled in while fitting).
            virtual const Core::Flux &mask_flux() const {return __mask_flux;}

            ///See __mask_flux (filled in while fitting).
            Core::Flux &mask_flux() {return __mask_flux;}

            ///The \f$\chi^2\f$ of the last scaling fit.
            virtual double &chi2() {return __chi2;}

            ///The reduced \f$\chi^2\f$ of the last scaling fit.
            double reduced_chi2() const
            {return __chi2/(flux_fit_pixel_count() - __group_sources);}

            ///Read-only \f$\chi^2\f$ of the last scaling fit.
            virtual double chi2() const {return __chi2;}

            ///The number of saturated pixels assigned to this source.
            virtual unsigned saturated_pixel_count() const =0;

            ///\brief the sum of (background excess)/(pixel variance+background
            ///variance) of the individual non-saturated pixels assigned to this
            ///source.
            virtual double merit() const = 0;

            ///The worst quality flag for any pixel assigned to the source.
            virtual Core::PhotometryFlag quality_flag() const = 0;

            ///The id the source was assigned in the source_assignment array.
            unsigned source_assignment_id() const
            {return __source_assignment_id;}

            ///Teh sub-pixel map supplied on input.
            const Core::SubPixelMap &subpix_map() const {return *__subpix_map;}

            ///Reject the source as non-point.
            void set_nonpoint() {__nonpoint=true;}

            ///Was the source previously rejected as non-point.
            bool is_nonpoint() const {return __nonpoint;}

            ///The reason the source was dropped, if it was.
            SourceDropReason drop_reason() const {return __drop_reason;}

            ///Drop this source for the given reason.
            void drop(SourceDropReason reason) {__drop_reason = reason;}

            ///The filename of the image from which this source was extracted.
            unsigned image_id() const {return __image_id;}

            ///The output filename where this source should be saved.
            const std::string &output_filename() const {return __output_filename;}

            ///The amplitude of the PSF
            double psf_amplitude() const {return __amplitude;}

            ///Sets the amplitude of the PSF to the given value.
            void set_psf_amplitude(double value) {__amplitude=value;}

            ///\brief Set how many sources were in the same group as this one
            ///when amplitude was fit.
            void set_sources_in_group(unsigned source_count)
            {__group_sources = source_count;}


            ///\brief Reference to the values of the expanion terms the PSF
            ///of the source depends on.
            Eigen::VectorXd &expansion_terms()
            {return __expansion_terms;}

            ///The values of the expanion terms the PSF of the source depends on.
            const Eigen::VectorXd &expansion_terms() const
            {return __expansion_terms;}

            ///\brief Fill a vector with the PSF parameters for this source given
            ///a set of polynomial expansion coefficients.
            void fill_psf_params(
                ///The polynomial expansion coefficients for all PSF terms.
                const Eigen::VectorXd &expansion_coef,

                ///The vector to fill with the PSF terms for this source. Resized
                ///as necessary.
                Eigen::VectorXd &psf_params
            );

        }; //End Source class.

    ///\brief Comparison between this and RHS, ordering by the
    ///respective merit function.
    template<class PSF_TYPE>
        bool operator<(const Source<PSF_TYPE> &lhs, const Source<PSF_TYPE> &rhs)
        {return lhs.merit()<rhs.merit();}

    ///Opposite of operator<()
    template<class PSF_TYPE>
        bool operator>(const Source<PSF_TYPE> &lhs, const Source<PSF_TYPE> &rhs)
        {return lhs.merit()>rhs.merit();}


    template<class PSF_TYPE>
        double Source<PSF_TYPE>::fit_flux(
            const Eigen::VectorXd &estimated_excesses
        )
        {
            double sum_sq_rhs,
                   sum_rhs_estimated;
            pixel_excess_reductions(estimated_excesses,
                                    sum_rhs_estimated,
                                    sum_sq_rhs);
            double sum_sq_estimated = estimated_excesses.squaredNorm(),
                   amplitude = sum_rhs_estimated / sum_sq_estimated,
                   result = flux(0).value() - amplitude;

            flux(0).value() = amplitude;
            chi2() = sum_sq_rhs - amplitude * sum_rhs_estimated;
            flux(0).error() = std::sqrt(chi2()
                                        /
                                        sum_sq_estimated
                                        /
                                        (flux_fit_pixel_count() - 1));

#ifdef VERBOSE_DEBUG
            std::cerr << "Individial flux fit for source ("
                      << this
                      << ") gave flux = "
                      << flux(0).value()
                      << std::endl;
#endif
            return result;
        }

    template<class PSF_TYPE>
        Source<PSF_TYPE>::Source(
            const Core::SourceID     &id,
            double                    x0,
            double                    y0,
            double                    gain,
            const Background::Source &background,
            int                       source_assignment_id,
            const Core::SubPixelMap  *subpix_map,
            unsigned                  image_id,
            const std::string        &output_fname
        )
        : PSF::MapSource(id, 1, x0, y0, background),
        __background_electrons(background.value() * gain),
        __background_electrons_variance(
            std::pow(background.error() * gain, 2)
        ),
        __gain(gain),
        __amplitude(Core::NaN),
        __source_assignment_id(source_assignment_id),
        __group_sources(1),
        __image_id(image_id),
        __subpix_map(subpix_map),
        __nonpoint(false),
        __mask_flux(Core::NaN),
        __output_filename(output_fname),
        __chi2(Core::NaN),
        __drop_reason(NOT_DROPPED)
    {}

    template<class PSF_TYPE>
        void Source<PSF_TYPE>::fill_psf_params(
            const Eigen::VectorXd &expansion_coef,
            Eigen::VectorXd &psf_params
        )
    {
        assert(expansion_coef.size() % __expansion_terms.size() == 0);

        unsigned num_poly_terms = __expansion_terms.size(),
                 num_psf_params = expansion_coef.size() / num_poly_terms;

        psf_params.resize(num_psf_params);
        for(
            unsigned psf_param_ind = 0;
            psf_param_ind < num_psf_params;
            ++psf_param_ind
        )
            psf_params[psf_param_ind] = __expansion_terms.dot(
                expansion_coef.segment(
                    psf_param_ind * num_poly_terms,
                    num_poly_terms
                )
            );
    }



} //End FitPSF namespace.

#endif
