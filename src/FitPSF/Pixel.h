/**\file
 *
 * \brief Defines a class describing a single pixel participating in PSF/PRF
 * fitting.
 *
 * \ingroup FitPSF
 */

#ifndef __PSF_FIT_PIXEL_H
#define __PSF_FIT_PIXEL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Background/Source.h"
#include "../Core/Typedefs.h"
#include <set>

namespace FitPSF {

    /**\brief A class tracking all pixel level informaiton for pixels
     * participating in PSF/PRF fitting.
     *
     * \ingroup FitPSF
     */
    template<class SOURCE_TYPE>
        class LIB_LOCAL Pixel {
        private:
            ///Alias type for a set of sources.
            typedef std::set< SOURCE_TYPE* > SourceSet;

            unsigned long
                ///The x-coordinate of the pixel within the image.
                __x,

                ///The y-coordinate of the pixel within the image.
                __y;

            double
                ///The measured value of the pixel in electrons.
                __measured,

                ///The estiamted variance of the pixel in electrons^2.
                __variance;

            ///The qualty flag for this pixel.
            Core::PhotometryFlag __flag;

            ///Is this pixel included in the PSF/PRF shape fitting.
            bool __shape_fit;

            ///\brief The index of this pixel within the flux fitting vector
            ///(negative if excluded).
            int __flux_fit_index;

            ///The set of sources this pixel belongs to.
            SourceSet __sources;

            ///Tell all sources that this pixel was excluded from the fit.
            void exclusion_notify_sources(
                ///Was the pixel previously included in the shape fit and is now
                ///excluded?
                bool from_shape_fit,

                ///Was the pixel previously included in the flux fit and is now
                ///excluded?
                bool from_flux_fit
            );

        public:
            ///Construct a new pixel belonging to a single source.
            Pixel(
                ///The x-coordinate of the pixel within the image.
                unsigned x,

                ///The y-coordinate of the pixel within the image.
                unsigned y,

                ///The measured value of the pixel and its variance in electrons
                ///and electrons^2 respectively.
                const std::pair<double, double> measured_variance,

                ///The quality flag for this pixel.
                Core::PhotometryFlag flag,

                ///The source this pixel belongs to. More sources can be added
                ///later.
                SOURCE_TYPE *source,

                ///The index within the flux fitting vector of this pixel (can be
                ///changed later).
                unsigned flux_fit_index = 0
            ) :
                __x(x),
                __y(y),
                __measured(measured_variance.first),
                __variance(measured_variance.second),
                __flag(flag),
                __shape_fit(true),
                __flux_fit_index(flux_fit_index)
            {assert(__variance > 0); add_to_source(source);}

            ///Add this pixel to another source.
            void add_to_source(
                ///The source to add the pixel to.
                SOURCE_TYPE *source
            )
            {
#ifdef VERBOSE_DEBUG
                std::cerr << "Adding source "
                          << source
                          << " to Pixel("
                          << x() << ", " << y()
                          << ") - has " << __sources.size()
                          << " sources:";
                for(
                    typename SourceSet::const_iterator
                        s_i = __sources.begin();
                    s_i != __sources.end();
                    ++s_i
                )
                    std::cerr << " " << *s_i;
                std::cerr << std::endl;
#endif
                __sources.insert(source);

#ifdef VERBOSE_DEBUG
                std::cerr << " Finished with "
                          << __sources.size()
                          << " sources:";
                for(
                    typename SourceSet::const_iterator
                        s_i = __sources.begin();
                    s_i != __sources.end();
                    ++s_i
                )
                    std::cerr << " " << *s_i;
                std::cerr << std::endl;
#endif
            }

            ///The set of sources this pixel belongs to.
            const SourceSet &sources() const
            {return __sources;}

            ///The x-coordinate of the pixel within the image.
            unsigned x() const {return __x;}

            ///The y-coordinate of the pixel within the image.
            unsigned y() const {return __y;}

            ///The measured value of the pixel in electrons.
            double measured() const {return __measured;}

            ///The estimated variance of the pixel in electrons^2.
            double variance() const {return __variance;}

            ///Should this pixel be included in the fit for the PSF shape?
            bool shape_fit() const {return __shape_fit;}

            ///\brief Should this pixel be included in fit for the amplitude of
            ///its sources?
            bool flux_fit() const {return __flux_fit_index >= 0;}

            ///Exclude this pixel from the fit for the PSF shape.
            void exclude_from_shape_fit()
            {
                if(__shape_fit) {
                    __shape_fit = false;
#ifdef VERBOSE_DEBUG
                    std::cerr << "Excluding pixel ("
                              << this
                              << ") at ("
                              << x()
                              << ", "
                              << y()
                              << ") from shape fitting."
                              << std::endl;
#endif
                    exclusion_notify_sources(true, false);
                }
            }

            ///Exclude this pixel from the fit for the PSF amplitudes.
            void exclude_from_flux_fit()
            {
                if(__flux_fit_index >= 0) {
                    __flux_fit_index = -1;
#ifdef VERBOSE_DEBUG
                    std::cerr << "Excluding pixel ("
                              << this
                              << ") at ("
                              << x()
                              << ", "
                              << y()
                              << ") from flux fitting."
                              << std::endl;
#endif
                    exclusion_notify_sources(false, true);
                }
            }

            ///\brief The index of this pixel within the flux fitting vector.
            ///
            ///Note that shape fitting pixels may have indenedent orbering from
            ///flux only fitting pixels.
            ///
            ///The behavior is undefined if the pixel is excluded from flux
            ///fitting.
            unsigned flux_fit_index() const
            {
                assert(__flux_fit_index >= 0);
                return static_cast<unsigned>(__flux_fit_index);
            }

            ///Set the index within the flux fitting vector for this pixel.
            void set_flux_fit_index(unsigned index)
            {assert(__flux_fit_index >= 0); __flux_fit_index = index;}

            ///Is the pixel saturated?
            Core::PhotometryFlag flag() const {return __flag;}

            ///Is the pixel shared by more than one source.
            bool shared() const {return __sources.size() > 1;}

        };//End Pixel class.

    ///\brief Return the square of the signal to noise ratio with which a
    ///pixel sticks above a background.
    double background_excess(
        ///The value of the pixel in electrons.
        double value,

        ///The variance of the pixel in electrons^2
        double variance,

        ///The background flux to assume in electrons.
        double background_value,

        ///The variance to assume for the background in electrons^2.
        double background_variance
    );

    ///See background_excess(double, double, dobule, double).
    template<class SOURCE_TYPE>
        double background_excess(
            ///The pixel for which to calculate the background excess.
            const Pixel<SOURCE_TYPE> &pixel,

            ///The background flux to assume in electrons.
            double background_value,

            ///The variance to assume for the background in electrons^2.
            double background_variance
        );

    ///See background_excess(double, double, double, double)
    double background_excess(
        ///The value of the pixel in electrons.
        double value,

        ///The variance of the pixel in electrons^2
        double variance,

        ///The background flux to assume in ADU.
        const Background::Source &background_adu,

        ///The gain to use when converting the background value from ADU to
        ///electrons.
        double gain
    );

    ///See background_excess(double, double, double, double).
    template<class SOURCE_TYPE>
        double background_excess(
            ///The pixel for which to calculate the background excess.
            const Pixel<SOURCE_TYPE> &pixel,

            ///The background flux in electrons to assume.
            const Background::Source &background,

            ///The gain to use when converting the background value from ADU to
            ///electrons.
            double gain
        );

    template<class SOURCE_TYPE>
        void Pixel<SOURCE_TYPE>::exclusion_notify_sources(bool from_shape_fit,
                                                          bool from_flux_fit)
        {
            for(
                typename SourceSet::iterator source_i = __sources.begin();
                source_i != __sources.end();
                ++source_i
            )
                (*source_i)->pixel_was_excluded(this,
                                                from_shape_fit,
                                                from_flux_fit);
        }

    template<class SOURCE_TYPE>
        double background_excess(const Pixel<SOURCE_TYPE>  &pixel,
                                 double                     background_value,
                                 double                     background_variance)
        {
            return background_excess(pixel.measured(),
                                     pixel.variance(),
                                     background_value,
                                     background_variance);
        }

    template<class SOURCE_TYPE>
        double background_excess(const Pixel<SOURCE_TYPE> &pixel,
                                 const Background::Source &background,
                                 double gain)
        {
            return background_excess(pixel.measured(),
                                     pixel.variance(),
                                     background,
                                     gain);
        }

} //End FitPSF namespace.

#endif
