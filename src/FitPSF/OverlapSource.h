/**\file
 *
 * \brief Defines the overlap-related features of PSF fitting sources.
 *
 * \ingroup FitPSF
 */

#ifndef __OVERLAP_SOURCE_H
#define __OVERLAP_SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Source.h"
#include "PixelOrder.h"
#include "Image.h"

namespace FitPSF {

    ///PSF fitting source which handles overlaps with other sources.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        class LIB_LOCAL OverlapSource : public Source<PSF_TYPE> {
        private:

            ///Alias type for a set of sources used in PSF fitting.
            typedef std::set< FIT_SOURCE_TYPE* > SourceSet;

            ///Alias type for a list of fit source pixels.
            typedef std::list< Pixel<FIT_SOURCE_TYPE>* > PixelList;

            ///Alias type for iterator to unmutable list of fit source pixels.
            typedef typename PixelList::const_iterator ConstPixelIter;

            ///Alias type for iterator to mutable list of fit source pixels.
            typedef typename PixelList::iterator PixelIter;

            ///\brief The pixels assigned to this source.
            ///
            ///Shape fitting pixels are first, followed by flux but not shape
            ///fitting and finally pixels excluded from both fits.
            PixelList __pixels;

            PixelIter
                ///\brief Iterator to the currently selected pixel (see
                ///restart_pixel_iteration() and next_pixel())
                __current_pixel,

                ///\brief Iterator to the first pixel that participates in flux
                ///fitting, but not shape fitting.
                __first_flux_fit_pixel,

                ///\brief Iterator to the first pixel that does not participate
                ///in either shape or flux fitting.
                __first_excluded_pixel;

            ConstPixelIter
                ///Unmutable version of __current_pixel.
                __const_current_pixel,

                ///Unmutable version of __first_flux_fit_pixel.
                __const_first_flux_fit_pixel,

                ///Unmutable version of __first_excluded_pixel.
                __const_first_excluded_pixel;

            ///\brief Some of the sources which this source overlaps with on the
            ///observed image.
            ///
            ///If a source is known to overlap, but none of the sources it
            ///overlaps with are known this set contains only the source's
            ///own id.
            SourceSet __overlaps;

            ///The worst flag for any of the pixels associated with the source.
            Core::PhotometryFlag __quality_flag;

            ///\brief The sum of the background excesses of all good
            ///source pixels not shared with another source.
            double __merit;

            ///\brief If non-zero the square of the largest allowed circular
            ///aperture for the source.
            double __max_circular_aperture2,

                   ///\brief The actual minimal circular aperture that contains
                   ///all source pixels squared.
                   __aperture2;

            unsigned
                ///The number of non-good pixels assigned to this source.
                __saturated_pixel_count,

                ///\brief The number of source pixels used in the last fit for
                ///the shape of the PSF.
                __shape_fit_pixel_count,

                ///\brief The number of pixels used in the last fit for the
                ///source flux.
                __flux_fit_pixel_count;

            ///Is this source ready for PSF fitting?
            bool __ready_to_fit;

            ///Assigns the given pixel to this source.
            void take_pixel(
                ///The x coordinate of the pixel to add
                unsigned x,

                ///The y coordinate of the pixel to add
                unsigned y,

                ///The PSF fitting image.
                Image<FIT_SOURCE_TYPE> &psffit_image
            )
            {
#ifdef VERBOSE_DEBUG
                std::cerr << "Source("
                          << this->x() << ", " << this->y() << "): "
                          << this
                          << " taking pixel(" << x << ", " << y << ")"
                          << std::endl;
#endif
                __pixels.push_back(
                    psffit_image.assign_to_source(
                        x,
                        y,
                        static_cast<FIT_SOURCE_TYPE*>(this)
                    )
                );
            }

            ///\brief Adds all pixels within the given aperture to the source
            ///taking care of overlaps.
            void add_pixels_in_circle(
                ///The size of the circle within which pixels are added to the
                ///source.
                double aperture2,

                ///The image used for the PSF fit.
                Image<FIT_SOURCE_TYPE> &psffit_image
            );

            ///\brief Adds all pixels within the given rectangle to the source
            ///taking care of overlaps.
            void add_pixels_in_rectangle(
                ///relative to source center
                double left,

                ///relative to source center
                double right,

                ///relative to source center
                double bottom,

                ///relative to source center
                double top,

                ///The Image used for fitting.
                Image<FIT_SOURCE_TYPE>& psffit_image
            );

            ///\brief Appropriately adds or not a pixel to a source.
            ///
            ///Checks if the given pixel is sufficiently above the background
            ///and if so checks for conflicts and if appropriate adds the
            ///pixel to this source, and calls itself on all neighboring
            ///pixels.
            double process_pixel(
                ///The x coordinate of the pixel to potentially add.
                unsigned               x,

                ///The y coordinate of the pixel to potentially add.
                unsigned               y,

                ///The PSF fitting image image.
                Image<FIT_SOURCE_TYPE> &psffit_image,

                ///See the alpha parameter of the constructor.
                double                 alpha,

                ///Is this the central pixel of the source.
                bool                   central_pixel = false
            );

            ///\brief Flag any pixels not suitable for fitting as excluded from
            ///shape/flux fits.
            void reject_pixels_from_fit();

            ///Return how constraining this source is for PSF fitting.
            void calculate_merit();

        protected:
            ///See Source::pixel_excess_reductions.
            void pixel_excess_reductions(
                const Eigen::VectorXd &to_dot_with,
                double &dot_product,
                double &excess_sum_squares
            );

        public:
            ///Creates a PSF fitting source by detecting bright pixels.
            OverlapSource(
                ///See same name argument to Source constructor.
                const Core::SourceID &id,

                ///See same name argument to Source constructor.
                double x0,

                ///See same name argument to Source constructor.
                double y0,

                ///See same name argument to Source constructor.
                const Background::Source &background,

                ///The image of fit pixels we are deriving the PSF/PRF map of. On
                ///exit, it is updated with the pixels belonging to the newly
                ///constructed source.
                Image<FIT_SOURCE_TYPE> &psffit_image,

                ///How much above the background a pixel needs to be in order to
                ///be allocated to this source (the alpha parameter in the
                ///description)
                double alpha,

                ///See same name argument to Source constructor.
                int source_assignment_id,

                ///See same name argument to Source constructor.
                const Core::SubPixelMap *subpix_map,

                ///If nonzero impose a circular aperture for the source no larger
                ///than the given value (otherwise uses only pixels inconsistent
                ///with the background at the prescribed by alpha level). The
                ///size of the circular aperture is the smallest size possible
                ///that encapsulates all pixels that pass the alpha test.
                double max_circular_aperture,

                ///See same name argument to Source constructor.
                const std::string &output_fname
            );

            ///Expose Source constructor (see for argument details).
            OverlapSource(
                const Core::SourceID &id,
                double x0,
                double y0,
                const Background::Source& background,
                Image<FIT_SOURCE_TYPE> &psffit_image,
                int source_assignment_id,
                const Core::SubPixelMap *subpix_map,
                double  left,
                double  right,
                double  bottom,
                double  top,
                const std::string &output_fname
            );

            ///\brief A post-processing step to be called after constructing
            ///all sources in the same image.
            virtual void finalize_pixels();

            ///See Source::ready_to_fit()
            bool ready_to_fit() const;

            ///\brief Adds the given set of overlaps to this sources overlaps,
            ///skipping this source if found.
            void add_overlaps(const SourceSet &extra_overlaps);

            ///\brief A set of the sources which this source overlaps with on the
            ///observed image.
            ///
            ///Until finalize_pixels() is called, this function returns an
            ///empty set!
            const SourceSet &overlaps() const {return __overlaps;}

            ///\brief Sets the entries in the flux fitting matrix
            ///corresponding to this source.
            template<class SHAPE_FIT_OUTPUT_TYPE, class FLUX_FIT_OUTPUT_TYPE>
                void fill_fluxfit_column(
                    ///A vector with the PSF parameters to assume.
                    const PSF::PSF &psf,

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

            ///\brief Sets the entries in the flux fitting matrix
            ///corresponding to this source assuming the given PSF map.
            template<class SHAPE_FIT_OUTPUT_TYPE, class FLUX_FIT_OUTPUT_TYPE>
                void fill_fluxfit_column(
                    ///A vector with the PSF parameters to assume.
                    const PSF::Map &psf_map,

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
                )
                {
                    PSF::PSF *psf = psf_map(*this);
                    fill_fluxfit_column(*psf,
                                        shape_fit_output,
                                        flux_fit_output,
                                        sequential_shape_fit_pixels,
                                        sequential_flux_fit_pixels);
                    delete psf;

                }

            using Source<PSF_TYPE>::fit_flux;

            ///\brief Sets the flux of the source to its PSF fitted value
            ///returning the change.
            ///
            ///This should only be used for isolated sources.
            ///
            ///\todo Handle saturated pixels.
            double fit_flux(
                ///The PSF whose amplitude we are fitting.
                const PSF::PSF &psf
            );

            ///\brief Return the integral of the PSF times sub-pixel map over
            ///the given pixel.
            virtual double pixel_psf(const Pixel<FIT_SOURCE_TYPE>* pixel,
                                     const PSF::PSF &psf) const;

            ///\brief The integral of the normalized PSF over the current pixel
            ///and its derivatives
            double pixel_psf(
                ///The PSF to integrate.
                const PSF::PSF &psf,

                ///The derivate order to return.
                PSF::SDKDerivative deriv = PSF::NO_DERIV
            ) const
            {
                if(deriv == PSF::NO_DERIV) return pixel_psf(*__current_pixel, psf);
                else throw Error::InvalidArgument("FitPSF::Source::pixel_psf",
                                                  "No derivatives allowed!");
            }

            ///\brief Return a const reference to the list of image pixels
            ///assigned to this source.
            const PixelList &pixels() const
            {return __pixels;}

            ///\brief Iterator over only shape fitting pixels, pointing to the
            ///first pixel included in shape fitting.
            PixelIter shape_fit_pixels_begin()
            {return __pixels.begin();}

            ///Constant version of shape_fit_pixels_begin().
            ConstPixelIter shape_fit_pixels_begin() const
            {return __pixels.begin();}

            ///\brief Iterator over only shape fitting pixels, pointing to one
            ///past the last pixel included in shape fitting.
            PixelIter shape_fit_pixels_end()
            {return __first_flux_fit_pixel;}

            ///\brief Constant version of shape_fit_pixels_end().
            ConstPixelIter shape_fit_pixels_end() const
            {return __const_first_flux_fit_pixel;}

            ///\brief Iterator over flux but not shape fitting pixels, pointing
            ///to the first pixel included in flux but not shape fitting.
            PixelIter flux_fit_pixels_begin()
            {return __first_flux_fit_pixel;}

            ///\brief Constant version of flux_fit_pixels_begin().
            ConstPixelIter flux_fit_pixels_begin() const
            {return __const_first_flux_fit_pixel;}

            ///\brief Iterator over flux but not shape fitting pixels, pointing
            ///to one past the last pixel included in flux but not shape fitting.
            PixelIter flux_fit_pixels_end()
            {return __first_excluded_pixel;}

            ///\brief Constant version of flux_fit_pixels_end().
            ConstPixelIter flux_fit_pixels_end() const
            {return __const_first_excluded_pixel;}

            ///The number of pixels assigned to this source
            size_t pixel_count() const {return __pixels.size();}

            ///See Source::saturated_pixel_count()
            unsigned saturated_pixel_count() const
            {return __saturated_pixel_count;}

            ///Restarts the iteration over pixels from the beginning
            virtual void restart_pixel_iteration()
            {__current_pixel = __pixels.begin();}

            ///\brief Advance to the next pixel of the source returning true
            ///if it is not past the last pixel.
            virtual bool next_pixel()
            {
                ++__current_pixel;
                ++__const_current_pixel;
                return __current_pixel != __pixels.end();
            }

            ///\brief The number of pixels belonging to this source suitable for
            ///PSF shape fitting.
            size_t shape_fit_pixel_count() const
            {return __shape_fit_pixel_count;}

            ///\brief The number of pixels belonging to this source suitable for
            ///PSF flux fitting.
            size_t flux_fit_pixel_count() const {return __flux_fit_pixel_count;}

            ///\brief Update the source to note that one of its pixels is now
            ///excluded from the fit.
            virtual void pixel_was_excluded(
                ///The pixel which was excluded. Must belong to this source,
                ///otherwise the behaviour is undefined. Calling the pixel's
                ///exclude_from_*_fit function(s) must be done before invoking
                ///this method.
                const Pixel<FIT_SOURCE_TYPE> *pixel,

                ///Was the pixel previously included in the shape fit and is now
                ///excluded?
                bool from_shape_fit,

                ///Was the pixel previously included in the flux fit and is now
                ///excluded?
                bool from_flux_fit
            );

            ///\brief Let the source know it was discarded from shape fitting.
            virtual void exclude_from_shape_fit();

            ///\brief Let the source know it was discarded from flux fitting.
            virtual void exclude_from_flux_fit();

            ///\brief Replace the source pixels.
            ///
            ///The iterators should point to elements that provide x(), y(),
            ///value(), variance() and saturated().
            template<class ITERATOR_TYPE>
                void replace_pixels(
                    ///Should point to the first pixel to insert.
                    ITERATOR_TYPE first,

                    ///Should point to one past the last pixel to insert.
                    ITERATOR_TYPE past_last
                );

            ///The quality flag for the current pixel.
            Core::PhotometryFlag current_pixel_flag() const
            {return (*__current_pixel)->flag();}

            ///The current pixel.
            const Pixel<FIT_SOURCE_TYPE> *current_pixel()
            {return *__const_current_pixel;}

            ///See Source::calculate_mask_flux()
            const Core::Flux &calculate_mask_flux(const PSF::PSF &psf);

            ///See Source::merit()
            double merit() const {assert(__ready_to_fit); return __merit;}

            ///The worst quality flag for any pixel assigned to the source.
            Core::PhotometryFlag quality_flag() const {return __quality_flag;}

            ///\brief The smallest circular aperture size that if centered on the
            ///source position contains all pixels assigned to the source.
            double aperture() const {return std::sqrt(__aperture2);}

        }; //End OverlapSource class.

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::add_pixels_in_circle(
            double aperture2,
            Image<FIT_SOURCE_TYPE> &psffit_image
        )
        {
            double yrange = std::sqrt(aperture2),
            dist_to_boundary = std::min(this->x(), this->y());
            dist_to_boundary = std::min(
                dist_to_boundary,
                std::min(
                    psffit_image.x_resolution() - this->x(),
                    psffit_image.y_resolution() - this->y()
                )
            );

            if(std::pow(dist_to_boundary, 2) < aperture2)
                __quality_flag = Core::BAD;
            long
                miny = std::max(
                    long(0),
                    static_cast<long>(std::ceil(this->y() - yrange - 0.5))
                ),
                     maxy = std::min(
                         static_cast<long>(psffit_image.y_resolution()
                                           -
                                           1),
                         static_cast<long>(std::floor(this->y() + yrange - 0.5))
                     );

            for(long y = miny; y <= maxy; y++) {
                double xrange = std::sqrt(aperture2
                                          -
                                          std::pow(0.5 + y-this->y(), 2));
                long
                    minx = std::max(
                        long(0),
                        static_cast<long>(std::ceil(this->x() - xrange - 0.5))
                    ),
                         maxx = std::min(
                             static_cast<long>(psffit_image.x_resolution()
                                               -
                                               1),
                             static_cast<long>(std::floor(this->x() + xrange - 0.5))
                         );

                for(long x = minx; x <= maxx; x++)
                    take_pixel(x, y, psffit_image);
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::add_pixels_in_rectangle(
            double left,
            double right,
            double bottom,
            double top,
            Image<FIT_SOURCE_TYPE>& psffit_image
        )
        {
            int xmin = std::floor( Source<PSF_TYPE>::x() + left );
            int xmax = std::ceil( Source<PSF_TYPE>::x() + right ) - 1;
            int ymin = std::floor( Source<PSF_TYPE>::y() + bottom );
            int ymax = std::ceil( Source<PSF_TYPE>::y() + top ) - 1;

            if( xmin < 0 ) {
                xmin = 0;
                __quality_flag = Core::BAD;
            }

            if( ymin < 0 ) {
                ymin = 0;
                __quality_flag = Core::BAD;
            }

            if(
                xmax >= static_cast<int>(
                    psffit_image.x_resolution()
                )
            ) {
                // last valid coord
                xmax = psffit_image.x_resolution() - 1;
                __quality_flag = Core::BAD;
            }

            if(
                ymax >= static_cast<int>(psffit_image.y_resolution())
            ) {
                // last valid coord
                ymax = psffit_image.y_resolution() - 1;
                __quality_flag = Core::BAD;
            }

#ifdef VERBOSE_DEBUG
            std::cerr << "Taking all pixels in range "
                << xmin << " <= x <= " << xmax
                << ", "
                << ymin << " <= y <= " << ymax
                << std::endl;
#endif
            for (int xx = xmin; xx <= xmax; ++xx) {
                for (int yy = ymin; yy <= ymax; ++yy) {
                    take_pixel(xx, yy, psffit_image);
                }
            }
#ifdef VERBOSE_DEBUG
            std::cerr << "Source ended up with " << pixel_count() << "pixels."
                << std::endl;
#endif
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        double OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::process_pixel(
            unsigned                    x,
            unsigned                    y,
            Image<FIT_SOURCE_TYPE>      &psffit_image,
            double                      alpha,
            bool                        central_pixel)
        {
            Core::PhotometryFlag
                pixel_flag = psffit_image.photometry_flag(x, y);
            __quality_flag = std::max(__quality_flag, pixel_flag);

            if(pixel_flag != Core::BAD || central_pixel) {
                double bg_excess = psffit_image.background_excess(
                    x,
                    y,
                    Source<PSF_TYPE>::background_electrons(),
                    Source<PSF_TYPE>::background_electrons_variance()
                );
                const Pixel<FIT_SOURCE_TYPE> *pixel = psffit_image(x, y);
                if(
                    !central_pixel
                    &&
                    (
                        (bg_excess < alpha)
                        ||
                        (
                            pixel != NULL
                            &&
                            (
                                pixel->shared()
                                ||
                                *(pixel->sources().begin()) != this
                            )
                        )
                    )
                )
                    return 0;
                else if(__max_circular_aperture2 == 0)
                    take_pixel(x, y, psffit_image);
            } else if(__max_circular_aperture2 == 0)
                take_pixel(x, y, psffit_image);

            unsigned long x_max = psffit_image.x_resolution() - 1,
                          y_max = psffit_image.y_resolution() - 1;
            double aperture2 = std::pow(0.5 + x - this->x(), 2)
                +
                std::pow(0.5 + y - this->y(), 2);
            if(x > 0)
                aperture2 = std::max(aperture2,
                                     process_pixel(x - 1,
                                                   y,
                                                   psffit_image,
                                                   alpha));
            else __quality_flag = Core::BAD;
            if(y>0) aperture2 = std::max(aperture2,
                                         process_pixel(x,
                                                       y - 1,
                                                       psffit_image,
                                                       alpha));
            else __quality_flag = Core::BAD;
            if(x < x_max) aperture2 = std::max(aperture2,
                                               process_pixel(x + 1,
                                                             y,
                                                             psffit_image,
                                                             alpha));
            else __quality_flag = Core::BAD;
            if(y < y_max) aperture2 = std::max(aperture2,
                                               process_pixel(x,
                                                             y + 1,
                                                             psffit_image,
                                                             alpha));
            else __quality_flag = Core::BAD;
            return aperture2;
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<
            FIT_SOURCE_TYPE,
            PSF_TYPE
        >::reject_pixels_from_fit()
        {
            for(
                PixelIter pix_i = __pixels.begin();
                pix_i != __pixels.end();
                ++pix_i
            ) {
                if((*pix_i)->flag() == Core::BAD) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Pixel("
                              << (*pix_i)->x()
                              << ", "
                              << (*pix_i)->y()
                              << ") is bad, excluding from shape and flux fit."
                              << std::endl;
#endif
                    (*pix_i)->exclude_from_shape_fit();
                    (*pix_i)->exclude_from_flux_fit();
                } else if(
                    (*pix_i)->shared()
                    ||
                    (*pix_i)->flag() != Core::GOOD
                    ||
                    std::isnan(this->background_electrons())
                    ||
                    std::isnan(this->background_electrons_variance())
                ) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Pixel("
                              << (*pix_i)->x()
                              << ", "
                              << (*pix_i)->y()
                              << ") is not perfect, excluding from shape fit."
                              << std::endl;
#endif
                    (*pix_i)->exclude_from_shape_fit();
                }
            }
            __pixels.sort(
                PixelOrder(Source<PSF_TYPE>::background_electrons(),
                           Source<PSF_TYPE>::background_electrons_variance())
            );
#ifdef VERBOSE_DEBUG
            std::cerr << "Ordered pixels:" << std::endl;
            for(
                PixelIter pix_i = __pixels.begin();
                pix_i != __pixels.end();
                ++pix_i
            ) {
                std::cerr << "Pixel("
                          << (*pix_i)->x()
                          << ", "
                          << (*pix_i)->y()
                          << ", shape: "
                          << (*pix_i)->shape_fit()
                          << ", flux: "
                          << (*pix_i)->flux_fit()
                          << ", merit: "
                          << background_excess(
                              **pix_i,
                              Source<PSF_TYPE>::background_electrons(),
                              Source<PSF_TYPE>::background_electrons_variance()
                          )
                          << ")"
                          << std::endl;
            }
#endif
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::calculate_merit()
        {
            __merit = 0.0;
            for(
                ConstPixelIter p = __pixels.begin();
                p != __pixels.end();
                ++p
            ) {
                if((*p)->flag() == Core::GOOD && !(*p)->shared()) {
                    __merit += background_excess(
                        **p,
                        Source<PSF_TYPE>::background_electrons(),
                        Source<PSF_TYPE>::background_electrons_variance()
                    );
                }
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::pixel_excess_reductions(
            const Eigen::VectorXd &to_dot_with,
            double &dot_product,
            double &excess_sum_squares
        )
        {
            dot_product = excess_sum_squares = 0;
            unsigned pix_ind = 0;
            for(
                ConstPixelIter p = __pixels.begin();
                p != __pixels.end();
                ++p
            ) {
                double bg_excess = background_excess(
                    **p,
                    Source<PSF_TYPE>::background_electrons(),
                    Source<PSF_TYPE>::background_electrons_variance()
                );
                dot_product += bg_excess * to_dot_with[pix_ind];
                excess_sum_squares += bg_excess * bg_excess;
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::OverlapSource(
            const Core::SourceID     &id,
            double                    x0,
            double                    y0,
            const Background::Source &background,
            Image<FIT_SOURCE_TYPE>   &psffit_image,
            double                    alpha,
            int                       source_assignment_id,
            const Core::SubPixelMap  *subpix_map,
            double                    max_circular_aperture,
            const std::string        &output_fname
        )
        : Source<PSF_TYPE>(id,
                           x0,
                           y0,
                           psffit_image.gain(),
                           background,
                           source_assignment_id,
                           subpix_map,
                           psffit_image.id(),
                           output_fname),
        __current_pixel(__pixels.begin()),
        __first_flux_fit_pixel(__pixels.end()),
        __first_excluded_pixel(__pixels.end()),
        __const_current_pixel(__pixels.begin()),
        __const_first_flux_fit_pixel(__pixels.end()),
        __const_first_excluded_pixel(__pixels.end()),
        __quality_flag(Core::GOOD),
        __merit(0),
        __max_circular_aperture2(max_circular_aperture > 0
                                 ? std::pow(max_circular_aperture, 2)
                                 : 0),
        __saturated_pixel_count(0),
        __shape_fit_pixel_count(0),
        __flux_fit_pixel_count(0),
        __ready_to_fit(false)
        {
            if(std::isnan(Source<PSF_TYPE>::background_electrons())) return;
            if(
                x0 < 0
                ||
                y0 < 0
                ||
                x0 >= psffit_image.x_resolution()
                ||
                y0 >= psffit_image.y_resolution()
            )
                return;

            __aperture2 = process_pixel(
                static_cast<unsigned long>(std::floor(x0)),
                static_cast<unsigned long>(std::floor(y0)),
                psffit_image,
                alpha,
                true
            );

            if(__max_circular_aperture2 != 0)
                add_pixels_in_circle(__aperture2, psffit_image);
            restart_pixel_iteration();
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::OverlapSource(
            const Core::SourceID      &id,
            double                     x0,
            double                     y0,
            const Background::Source  &background,
            Image<FIT_SOURCE_TYPE>    &psffit_image,
            int                        source_assignment_id,
            const Core::SubPixelMap   *subpix_map,
            double                     left,
            double                     right,
            double                     bottom,
            double                     top,
            const std::string         &output_fname
        )
        : Source<PSF_TYPE>(id,
                           x0,
                           y0,
                           psffit_image.gain(),
                           background,
                           source_assignment_id,
                           subpix_map,
                           psffit_image.id(),
                           output_fname),
        __current_pixel(__pixels.begin()),
        __first_flux_fit_pixel(__pixels.end()),
        __first_excluded_pixel(__pixels.end()),
        __const_current_pixel(__pixels.begin()),
        __const_first_flux_fit_pixel(__pixels.end()),
        __const_first_excluded_pixel(__pixels.end()),
        __quality_flag(Core::GOOD),
        __merit(0),
        __saturated_pixel_count(0),
        __shape_fit_pixel_count(0),
        __flux_fit_pixel_count(0),
        __ready_to_fit(false)
    {
        if (
            x0 < 0
            ||
            y0 < 0
            ||
            x0 >= psffit_image.x_resolution()
            ||
            y0 >= psffit_image.y_resolution()
        ) {
            // source is outside boundaries
            return;
        }

        add_pixels_in_rectangle(left,
                                right,
                                bottom,
                                top,
                                psffit_image);
        restart_pixel_iteration();
    }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::finalize_pixels()
        {
            if(__ready_to_fit) return;
#ifdef VERBOSE_DEBUG
            std::cerr << "Finalizing pixels of source ("
                      << this
                      << ") at ("
                      << this->x()
                      << ", "
                      << this->y()
                      << ")"
                      << std::endl;
#endif
            reject_pixels_from_fit();

            __current_pixel = __pixels.begin();

            __shape_fit_pixel_count =
                __flux_fit_pixel_count =
                __saturated_pixel_count = 0;
            for(
                __first_flux_fit_pixel = __current_pixel;
                (
                    __first_flux_fit_pixel != __pixels.end()
                    &&
                    (*__first_flux_fit_pixel)->shape_fit()
                );
                ++__first_flux_fit_pixel
            ) {
                (*__first_flux_fit_pixel)->set_flux_fit_index(
                    __shape_fit_pixel_count
                );
                ++__shape_fit_pixel_count;
                ++__flux_fit_pixel_count;
                if((*__first_flux_fit_pixel)->flag() == Core::SATURATED)
                    ++__saturated_pixel_count;

#ifdef VERBOSE_DEBUG
                std::cerr
                    << "Pixel ("
                    << *__first_flux_fit_pixel
                    << ") at ("
                    << (*__first_flux_fit_pixel)->x()
                    << ", "
                    << (*__first_flux_fit_pixel)->y()
                    << ") is "
                    << ((*__first_flux_fit_pixel)->shared() ? "" : "not ")
                    << "shared and will "
                    << ((*__first_flux_fit_pixel)->shape_fit() ? "" : "not ")
                    << "participate in shape fitting."
                    << std::endl;
#endif
                if((*__first_flux_fit_pixel)->shared())
                    add_overlaps((*__first_flux_fit_pixel)->sources());
            }

            for(
                __first_excluded_pixel = __first_flux_fit_pixel;
                (
                    __first_excluded_pixel != __pixels.end()
                    &&
                    (*__first_excluded_pixel)->flux_fit()
                );
                ++__first_excluded_pixel
            ) {
                assert(!(*__first_excluded_pixel)->shape_fit());
                if(!(*__first_excluded_pixel)->shared()) {
                    (*__first_excluded_pixel)->set_flux_fit_index(
                        __flux_fit_pixel_count
                        -
                        __shape_fit_pixel_count
                    );
#ifdef VERBOSE_DEBUG
                    std::cerr << "Source pixel("
                        << (*__first_excluded_pixel)->x()
                        << ", "
                        << (*__first_excluded_pixel)->y()
                        << ") assigned to flux fit index "
                        << (*__first_excluded_pixel)->flux_fit_index()
                        << std::endl;
#endif
                }
                if((*__first_excluded_pixel)->flag() == Core::SATURATED)
                    ++__saturated_pixel_count;

                if((*__first_excluded_pixel)->shared())
                    add_overlaps((*__first_excluded_pixel)->sources());

                ++__flux_fit_pixel_count;
            }
            assert(__first_excluded_pixel == __pixels.end()
                   ||
                   !(*__first_excluded_pixel)->flux_fit());
            __const_current_pixel = __current_pixel;
            __const_first_flux_fit_pixel = __first_flux_fit_pixel;
            __const_first_excluded_pixel = __first_excluded_pixel;
            calculate_merit();
            __ready_to_fit = true;
        }

    ///Is this source ready for PSF fitting?
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        bool OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::ready_to_fit() const
        {
            return (
                __ready_to_fit
                &&
                Source< PSF::PiecewiseBicubic >::ready_to_fit()
            );
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::add_overlaps(
            const SourceSet &extra_overlaps
        )
        {
            for(
                typename SourceSet::const_iterator
                    s = extra_overlaps.begin();
                s != extra_overlaps.end();
                ++s
            )
                if(*s != this) __overlaps.insert(*s);
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        template<class SHAPE_FIT_OUTPUT_TYPE, class FLUX_FIT_OUTPUT_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::fill_fluxfit_column(
            const PSF::PSF &psf,
            SHAPE_FIT_OUTPUT_TYPE shape_fit_output,
            FLUX_FIT_OUTPUT_TYPE flux_fit_output,
            bool sequential_shape_fit_pixels,
            bool sequential_flux_fit_pixels
        )
        {
            assert(shape_fit_output.size() == __shape_fit_pixel_count);
            assert(flux_fit_output.size() == (__flux_fit_pixel_count
                                              -
                                              __shape_fit_pixel_count));
            unsigned pix_ind = 0;
            for(
                ConstPixelIter pix_i = __pixels.begin();
                pix_i != __first_excluded_pixel;
                ++pix_i
            ) {
                double excess = (
                    pixel_psf(*pix_i, psf)
                    /
                    std::sqrt((*pix_i)->variance()
                              +
                              Source<PSF_TYPE>::background_electrons_variance())
                );
                if(pix_ind < __shape_fit_pixel_count)
                    shape_fit_output(
                        (
                            sequential_shape_fit_pixels
                            ? pix_ind
                            : (*pix_i)->flux_fit_index()
                        ),
                        1
                    ) = excess;
                else
                    flux_fit_output(
                        (
                            sequential_flux_fit_pixels
                            ? pix_ind - __shape_fit_pixel_count
                            : (*pix_i)->flux_fit_index()
                        ),
                        1
                    ) = excess;
                ++pix_ind;
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        double OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::fit_flux(
            const PSF::PSF &psf
        )
        {
            Eigen::VectorXd estimated_excesses(__flux_fit_pixel_count);
            Eigen::VectorBlock<Eigen::VectorXd>
                shape_fit_excesses = estimated_excesses.head(
                    __shape_fit_pixel_count
                ),
                flux_fit_excesses = estimated_excesses.tail(
                    __flux_fit_pixel_count
                    -
                    __shape_fit_pixel_count
                );
            fill_fluxfit_column(
                psf,
                shape_fit_excesses,
                flux_fit_excesses
            );
            return fit_flux(estimated_excesses);
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        double OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::pixel_psf(
            const Pixel<FIT_SOURCE_TYPE>* pixel,
            const PSF::PSF &psf
        ) const
        {
            double y_step = (
                1.0
                /
                std::max(
                    Source<PSF_TYPE>::subpix_map().y_resolution(),
                    ulong1
                )
            );
            double x_step = (
                1.0
                /
                std::max(
                    Source<PSF_TYPE>::subpix_map().x_resolution(),
                    ulong1
                )
            );

            assert(pixel->flux_fit());
            double y0 = pixel->y() + 0.5 * y_step - this->y();
            if(
                Source<PSF_TYPE>::subpix_map().x_resolution() == 0
                &&
                Source<PSF_TYPE>::subpix_map().y_resolution() == 0
            ) {
                return psf(pixel->x() + 0.5 * x_step - this->x(), y0);
            } else {
                double result = 0.0;
                for(
                    unsigned subpix_y = 0;
                    subpix_y < Source<PSF_TYPE>::subpix_map().y_resolution();
                    ++subpix_y
                ) {
                    double x0 = pixel->x() + 0.5 * x_step - this->x();
                    for(
                        unsigned subpix_x = 0;
                        subpix_x
                        <
                        Source<PSF_TYPE>::subpix_map().x_resolution();
                        ++subpix_x
                    ) {
                        result += (
                            Source<PSF_TYPE>::subpix_map()(subpix_x, subpix_y)
                            *
                            psf.integrate(x0, y0, x_step, y_step)
                        );
                        x0 += x_step;
                    }
                    y0 += y_step;
                }
                return result;
            }
        }


    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::pixel_was_excluded(
            const Pixel<FIT_SOURCE_TYPE> *pixel,
            bool from_shape_fit,
            bool from_flux_fit
        )
        {
            if(!__ready_to_fit) return;
#ifndef NDEBUG
            std::cerr << "Source at (x=" << this->x() << ", y=" << this->y()
                      << ") unreadied to fit" << std::endl;
#endif
            if(from_shape_fit) {
                assert(!pixel->shape_fit());
                __merit -= background_excess(
                    *pixel,
                    Source<PSF_TYPE>::background_electrons(),
                    Source<PSF_TYPE>::background_electrons_variance()
                );
                --__shape_fit_pixel_count;
                __ready_to_fit = false;
            }
            if(from_flux_fit) {
                assert(!pixel->flux_fit());
                --__flux_fit_pixel_count;
                __ready_to_fit = false;
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::exclude_from_shape_fit()
        {
#ifdef VERBOSE_DEBUG
            std::cerr << "Excluding source ("
                      << this
                      << ") at ("
                      << this->x()
                      << ", "
                      << this->y()
                      << ") from shape fit."
                      << std::endl;
#endif
            for(
                PixelIter pix_i = shape_fit_pixels_begin();
                pix_i != shape_fit_pixels_end();
                ++pix_i
            ) {
                assert(!(*pix_i)->shared());
                (*pix_i)->exclude_from_shape_fit();
                pixel_was_excluded(*pix_i, true, false);
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::exclude_from_flux_fit()
        {
#ifdef VERBOSE_DEBUG
            std::cerr << "Excluding source ("
                      << this
                      << ") at ("
                      << this->x()
                      << ", "
                      << this->y()
                      << ") from flux fit."
                      << std::endl;
#endif
            for(
                PixelIter pix_i = __pixels.begin();
                pix_i != __pixels.end();
                ++pix_i
            ) {
                (*pix_i)->exclude_from_flux_fit();
                pixel_was_excluded(*pix_i, false, true);
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        template<class ITERATOR_TYPE>
        void OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::replace_pixels(
            ITERATOR_TYPE first,
            ITERATOR_TYPE past_last
        )
        {
            __pixels.clear();
            __saturated_pixel_count = 0;
            __shape_fit_pixel_count = 0;
            __pixels.insert(__pixels.end(), first, past_last);
            for( ; first != past_last; ++first) {
                if(first->saturated()) ++__saturated_pixel_count;
                else if(!first->shared()) ++__shape_fit_pixel_count;
            }
        }

    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        const Core::Flux &OverlapSource<FIT_SOURCE_TYPE, PSF_TYPE>::calculate_mask_flux(
            const PSF::PSF &psf
        )
        {
            double captured_flux_fraction = 0,
                   measured_flux = 0,
                   measured_variance = 0;
            for(
                ConstPixelIter pix_i = __pixels.begin();
                pix_i != __pixels.end();
                ++pix_i
            ) {
                if((*pix_i)->flag() != Core::GOOD) continue;
                double pixel_left = (*pix_i)->x() - Source<PSF_TYPE>::x();
                double pixel_bottom = (*pix_i)->y() - Source<PSF_TYPE>::y();
                if(
                    Source<PSF_TYPE>::subpix_map().x_resolution() == 0
                    &&
                    Source<PSF_TYPE>::subpix_map().y_resolution() == 0
                ) {
                    captured_flux_fraction += psf(pixel_left + 0.5,
                                                  pixel_bottom + 0.5);
                } else {
                    double
                        x_step = (
                            1.0
                            /
                            Source<PSF_TYPE>::subpix_map().x_resolution()
                        );
                    double y_step = (
                        1.0
                        /
                        Source<PSF_TYPE>::subpix_map().y_resolution()
                    );
                    double y0 = pixel_bottom + 0.5 * y_step;
                    for(
                        unsigned subpix_y = 0;
                        subpix_y
                        <
                        Source<PSF_TYPE>::subpix_map().y_resolution();
                        ++subpix_y
                    ) {
                        double x0 = pixel_left + 0.5 * x_step;
                        for(
                            unsigned subpix_x = 0;
                            subpix_x
                            <
                            Source<PSF_TYPE>::subpix_map().x_resolution();
                            ++subpix_x
                        ) {
                            captured_flux_fraction += (
                                Source<PSF_TYPE>::subpix_map()(subpix_x, subpix_y)
                                *
                                psf.integrate(x0, y0, x_step, y_step)
                            );
                            x0 += x_step;
                        }
                        y0 += y_step;
                    }
                }
                measured_flux += (*pix_i)->measured();
                measured_variance += (*pix_i)->variance();
            }
            Source<PSF_TYPE>::mask_flux().value() = (measured_flux
                                                     /
                                                     captured_flux_fraction);
            Source<PSF_TYPE>::mask_flux().error() = (
                std::sqrt(measured_variance)
                /
                captured_flux_fraction
            );
            return Source<PSF_TYPE>::mask_flux();
        }



} //End FitPSF namespace.

#endif
