/**\file
 *
 * \brief Declares a class used for PSF fitting where pixel responses are
 * linear functions of the shape parameters.
 *
 * \ingroup FitPSF
 */

#ifndef __LINEAR_SOURCE_H
#define __LINEAR_SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "OverlapSource.h"
#include "../PSF/PiecewiseBicubic.h"
#include "../Background/Source.h"
#include "../Core/SourceID.h"
#include "../Core/Flux.h"
#include "Eigen/Dense"
#include <vector>
#include <valarray>
#include <list>

namespace FitPSF {

    ///\brief A class for PSF fitting sources whose flux distribution depends
    ///linearly on the PSF paramaters (for a fixed amplitude).
    class LIB_PUBLIC LinearSource :
        public OverlapSource<LinearSource, PSF::PiecewiseBicubic> {
    private:

        ///Convenience alias for unmutable iterators to source pixels.
        typedef std::list< Pixel<LinearSource>* >::const_iterator
            ConstPixelIter;

        ///Convenience alias for mutable iterators to source pixels.
        typedef std::list< Pixel<LinearSource>* >::const_iterator
            PixelIter;

        ///See psf argument of the constructors.
        const PSF::PiecewiseBicubic &__psf;

        ///\brief A (piece of a) matrix which after applying to the PSF
        ///fitting coefficients and scaling by the amplitude gives the
        ///integral of the PSF over the source pixels selected for shape
        ///fitting.
        const Eigen::Block<Eigen::MatrixXd> *__shape_fit_integral_matrix;

        ///\brief A (piece of a) vector which after scaling by the amplitude
        ///gives the integral of a PSF with only the overall integral
        ///parameters for the source pixels selected for shape fitting.
        const Eigen::VectorBlock<Eigen::VectorXd> *__shape_fit_offset;

        ///\brief A (piece of a) vector of the background excess values of
        ///the source pixels used in shape fitting.
        const Eigen::VectorBlock<Eigen::VectorXd> *__shape_fit_bg_excess;

        ///\brief The same as __shape_fit_integral_matrix, but for pixels
        ///participating in flux fitting only (excluded from shape fitting).
        Eigen::MatrixXd __flux_fit_integral_matrix;

        ///\brief The same as __shape_fit_offset, but for pixels
        ///participating in flux fitting only (excluded from shape fitting).
        Eigen::VectorXd __flux_fit_offset;

        ///\brief The same as __shape_fit_bg_excess, but for pixels
        ///participating in flux fitting only (excluded from shape fitting).
        Eigen::VectorXd __flux_fit_bg_excess;

        ///Is this source ready for PSF fitting.
        bool __ready_to_fit;

        ///\brief Calculate the integral of PSF * sub-pixel map over the
        ///given pixel for a set of PSF paramaters.
        void calculate_predicted_pixel_values(
            ///The left boundary of the pixel
            double pixel_left,

            ///The bottom bounary of the pixel
            double pixel_bottom,

            ///The collection of PSF parameters for which to calculate the
            ///predicted pixel response.
            const std::vector< std::valarray<double> > &parameter_sets,

            ///The array to fill with the final result. Must already have the
            ///correct size.
            std::valarray<double> &pixel_integrals
        );

        ///\brief Prepare matrices giving PSF integrals over the source
        ///pixels.
        ///
        ///The result is split between two matrices: the output argument
        ///to this method, with entries for shape fitting pixels only and
        ///an internally stored matrix for flux-only fitting pixels.
        void fill_pixel_integral_matrix(
            ///See same name argument of prepare_for_fitting.
            const std::vector< std::valarray<double> > &basis_parameter_sets,

            ///See same name argument of prepare_for_fitting.
            Eigen::Block<Eigen::MatrixXd> &shape_fit_integral_matrix,

            ///See same name argument of prepare_for_fitting.
            Eigen::VectorBlock<Eigen::VectorXd> &shape_fit_offset
        );

        ///\brief Fill vectors with the background excess values of the
        ///source pixels used for fitting.
        ///
        ///Just as fill_pixel_integral_matrix(), the result is split into two
        ///pieces: the output argument to this method, with entries for shape
        ///fitting pixels only and an internally stored matrix for flux-only
        ///fitting pixels.
        void fill_background_excess(
            ///See same name argument of prepare_for_fitting.
            Eigen::VectorBlock<Eigen::VectorXd> &shape_fit_background_excess
        );

        ///\brief Re-orders a flux fitting column from the order matching the
        ///list of source pixels to the order specified by the flux_fit_index()
        ///Pixel method.
        template<class OUTPUT_TYPE>
            void reorder_flux_fit_column(
                ///The input vector to re-order. The entries should be in the
                ///same order as the sub-set of source pixels (shape vs flux fit)
                ///for which this vector applies.
                const Eigen::VectorXd &input,

                ///Are the pixels being re-ordered shape fitting (vs flux
                ///fitting only)?
                bool shape_fit,

                ///The output vector to fill. Entries are according to the
                ///flux_fit_index() of the pixels.
                OUTPUT_TYPE &output
            );

    protected:
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
        );

    public:
        ///\brief See Source constructor with matching arguments (except the
        ///first two).
        LinearSource(
            ///A properly constructed PSF for this source. The exact shape
            ///parameters set are irrelevant. Must not be destroyed while
            ///this object is in use.
            const PSF::PiecewiseBicubic &psf,

            ///See OverlapSource::OverlapSource()
            const Core::SourceID &id,

            ///See OverlapSource::OverlapSource()
            double x0,

            ///See OverlapSource::OverlapSource()
            double y0,

            ///See OverlapSource::OverlapSource()
            const Background::Source &background,

            ///See OverlapSource::OverlapSource()
            Image<LinearSource> &psffit_image,

            ///See OverlapSource::OverlapSource()
            double alpha,

            ///See OverlapSource::OverlapSource()
            int source_assignment_id,

            ///See OverlapSource::OverlapSource()
            const Core::SubPixelMap *subpix_map,

            ///See OverlapSource::OverlapSource()
            double max_circular_aperture,

            ///See OverlapSource::OverlapSource()
            const std::string &output_fname
        ) :
            OverlapSource<LinearSource, PSF::PiecewiseBicubic>(
                id,
                x0,
                y0,
                background,
                psffit_image,
                alpha,
                source_assignment_id,
                subpix_map,
                max_circular_aperture,
                output_fname
            ),
            __psf(psf),
            __shape_fit_integral_matrix(NULL),
            __shape_fit_offset(NULL),
            __shape_fit_bg_excess(NULL),
            __ready_to_fit(false)
        {
#ifdef VERBOSE_DEBUG
            std::cerr << "Created source at " << this << std::endl;
#endif
        }

        ///\brief See Source constructor with matching arguments (except the
        ///first one).
        LinearSource(
            ///A properly constructed PSF for this source. The exact shape
            ///parameters set are irrelevant. Must not be destroyed while
            ///this object is in use.
            const PSF::PiecewiseBicubic &psf,

            ///See OverlapSource::OverlapSource()
            const Core::SourceID &id,

            ///See OverlapSource::OverlapSource()
            double x0,

            ///See OverlapSource::OverlapSource()
            double y0,

            ///See OverlapSource::OverlapSource()
            const Background::Source& background,

            ///See OverlapSource::OverlapSource()
            Image<LinearSource> &psffit_image,

            ///See OverlapSource::OverlapSource()
            int source_assignment_id,

            ///See OverlapSource::OverlapSource()
            const Core::SubPixelMap *subpix_map,

            ///See OverlapSource::OverlapSource()
            double left,

            ///See OverlapSource::OverlapSource()
            double right,

            ///See OverlapSource::OverlapSource()
            double bottom,

            ///See OverlapSource::OverlapSource()
            double top,

            ///See OverlapSource::OverlapSource()
            const std::string &output_fname
        ) :
            OverlapSource<LinearSource, PSF::PiecewiseBicubic>(
                id,
                x0,
                y0,
                background,
                psffit_image,
                source_assignment_id,
                subpix_map,
                left,
                right,
                bottom,
                top,
                output_fname
            ),
            __psf(psf),
            __shape_fit_integral_matrix(NULL),
            __shape_fit_offset(NULL),
            __shape_fit_bg_excess(NULL),
            __ready_to_fit(false)
        {
#ifdef VERBOSE_DEBUG
            std::cerr << "Created source at " << this << std::endl;
#endif
        }

        ///\brief A post-processing step to prepare the source for fitting,
        ///which must be called after all sources are constructed.
        ///
        ///In addition to setting up the source it also fills pieces of the
        ///matrices needed for fitting.
        virtual void prepare_for_fitting(
            ///The collection of PSF parameter basis vectors for which to
            ///calculate the predicted pixel response. The first set should
            ///correspond to the overall integral being projected out of the
            ///basis.
            const std::vector< std::valarray<double> > &basis_parameter_sets,

            ///A matrix piece filled by this method giving the integrals
            ///over the shape fitting source pixels but with only one
            ///basis vector contributing (\f$\mathbf{\tilde{M}^i}\f$ in
            ///the documentation). Must already have the correct size.
            ///The order of the rows is the same as the order of source
            ///pixels, but with non shape-fitting pixels skipped.
            Eigen::Block<Eigen::MatrixXd> shape_fit_integral_matrix,

            ///A vector piece filled by this method, giving the offset to
            ///apply during the fitting due to the basis vectors selected to
            ///have zero overall PSF integral.
            Eigen::VectorBlock<Eigen::VectorXd> shape_fit_offset,

            ///The vector to fill with the background excess values for shape
            ///fitting pixels.
            Eigen::VectorBlock<Eigen::VectorXd> shape_fit_background_excess
        );

        ///Is this source ready for PSF fitting?
        bool ready_to_fit() const
        {
            return (
                __ready_to_fit
                &&
                OverlapSource<
                    LinearSource,
                    PSF::PiecewiseBicubic
                >::ready_to_fit()
            );
        }

        using OverlapSource<LinearSource,
                            PSF::PiecewiseBicubic>::fill_fluxfit_column;

        ///\brief Sets the entries in the flux fitting matrix
        ///corresponding to this source.
        template<class SHAPE_FIT_OUTPUT_TYPE, class FLUX_FIT_OUTPUT_TYPE>
            void fill_fluxfit_column(
                ///A vector with the PSF expansion coefficients to assume.
                const Eigen::VectorXd &psf_expansion_coef,

                ///The location to fill with the predicted pixel responses
                ///for the given PSF parameters for shape fitting pixels.
                ///The order follows the order of source pixels. Must already
                ///have the correct size.
                SHAPE_FIT_OUTPUT_TYPE shape_fit_output,

                ///The location to fill with the predicted pixel responses
                ///for the given PSF parameters for flux but not shape
                ///fitting pixels. Which entry gets filled is determined by
                ///each pixel's flux_fit_index, so pass the entire vector
                ///rather than just the segment for this source. Must already
                ///have the correct size.
                FLUX_FIT_OUTPUT_TYPE flux_fit_output,

                ///Is it safe to assume that shape fitting pixels are
                ///sequentially ordered in shape_fit_output?
                bool sequential_shape_fit_pixels = true,

                ///Is it safe to assume that flux fitting only pixels are
                ///sequentially ordered in shape_fit_output?
                bool sequential_flux_fit_pixels = true

            );

        using OverlapSource<LinearSource, PSF::PiecewiseBicubic>::pixel_psf;

        ///\brief The integral of the normalized PSF over the current
        ///pixel and its derivatives
        double pixel_psf(PSF::SDKDerivative =PSF::NO_DERIV) const
        {
            throw Error::Fitting("PSF integrals should never be "
                                 "requested for piecewise PSF fits.");
        }

        using OverlapSource<LinearSource, PSF::PiecewiseBicubic>::fit_flux;

        ///\brief Fit for the flux of the source assuming the given PSF shape
        ///parameters and return the change from its previous value.
        double fit_flux(
            ///A vector with the expansion coefficients of the PSF parameters
            ///to assume.
            const Eigen::VectorXd &psf_expansion_coef
        );

        ///The signal to noise ratio.
        double signal_to_noise() const {return std::sqrt(merit());}

        ///Clean-up allocated fitting matrices.
        ~LinearSource()
        {
            if(__shape_fit_integral_matrix) delete __shape_fit_integral_matrix;
            if(__shape_fit_offset) delete __shape_fit_offset;
            if(__shape_fit_bg_excess) delete __shape_fit_bg_excess;
        }

    }; //End LinearSource class.

    template<class OUTPUT_TYPE>
        void LinearSource::reorder_flux_fit_column(
            const Eigen::VectorXd &input,
            bool shape_fit,
            OUTPUT_TYPE &output
        )
        {
            ConstPixelIter first_pix, last_pix;
            if(shape_fit) {
                first_pix = shape_fit_pixels_begin();
                last_pix = shape_fit_pixels_end();
            } else {
                first_pix = flux_fit_pixels_begin();
                last_pix = flux_fit_pixels_end();
            }

            unsigned input_index = 0;
            for(
                ConstPixelIter pix_i = first_pix;
                pix_i != last_pix;
                ++pix_i
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Setting flux fit index "
                          << (*pix_i)->flux_fit_index() << std::endl;
#endif
                output((*pix_i)->flux_fit_index(), 0) = input[input_index++];
            }
            assert(input_index == input.size());
        }

    template<class SHAPE_FIT_OUTPUT_TYPE, class FLUX_FIT_OUTPUT_TYPE>
        void LinearSource::fill_fluxfit_column(
            const Eigen::VectorXd &psf_expansion_coef,
            SHAPE_FIT_OUTPUT_TYPE shape_fit_output,
            FLUX_FIT_OUTPUT_TYPE flux_fit_output,
            bool sequential_shape_fit_pixels,
            bool sequential_flux_fit_pixels
        )
        {
            assert(ready_to_fit());

            assert(shape_fit_output.size()
                   ==
                   static_cast<int>(shape_fit_pixel_count()));
            assert(flux_fit_output.size()
                   >=
                   static_cast<int>(flux_fit_pixel_count()
                                    -
                                    shape_fit_pixel_count()));

            Eigen::VectorXd psf_params;
            fill_psf_params(psf_expansion_coef, psf_params);

#ifdef VERBOSE_DEBUG
            std::cerr << "Shape fit matrix ("
                << __shape_fit_integral_matrix->rows()
                << "x"
                << __shape_fit_integral_matrix->cols()
                << "): " << std::endl
                << *__shape_fit_integral_matrix
                << std::endl;

            std::cerr << "Flux fit matrix ("
                << __flux_fit_integral_matrix.rows()
                << "x"
                << __flux_fit_integral_matrix.cols()
                << "): " << std::endl
                << __flux_fit_integral_matrix
                << std::endl;
            std::cerr << "PSF params: " << psf_params << std::endl;
#endif

            if(sequential_shape_fit_pixels) {
                shape_fit_output = (
                    (*__shape_fit_integral_matrix) * psf_params
                    +
                    (*__shape_fit_offset)
                );
            } else {
                Eigen::VectorXd scrambled = (
                    (*__shape_fit_integral_matrix) * psf_params
                    +
                    (*__shape_fit_offset)
                );
                reorder_flux_fit_column(scrambled, true, shape_fit_output);
            }
            if(sequential_flux_fit_pixels) {
                flux_fit_output = (
                    __flux_fit_integral_matrix * psf_params
                    +
                    __flux_fit_offset
                );
            } else {
                Eigen::VectorXd scrambled = (
                    __flux_fit_integral_matrix * psf_params
                    +
                    __flux_fit_offset
                );
                reorder_flux_fit_column(scrambled, false, flux_fit_output);
            }
        }

} //End FitPSF namespace.

#endif
