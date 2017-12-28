/**\file
 *
 * \brief Defines some of the functions needed exclusively by the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#ifndef __PSF_FITTING_H
#define __PSF_FITTING_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Background/Source.h"
#include "../Background/Measure.h"
#include "../Core/SourceLocation.h"
#include "../Core/Image.h"
#include "../Core/SubPixelMap.h"
#include "../PSF/PiecewiseBicubic.h"
#include "Eigen/Dense"
#include <list>
#include <set>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_poly.h>
#include <iostream>

namespace FitPSF {

    template <class FIT_SOURCE_TYPE> class Image;
    template <class PSF_TYPE> class Source;

    class LinearSource;

    typedef std::list<LinearSource *> LinearSourceList;

    const unsigned long ulong1 = 1;

    /**\brief Define tags for reasons to exclude sources from the fit.
     *
     * \ingroup FitPSF
     */
    enum LIB_PUBLIC SourceDropReason {
        FEW_PIXELS,       ///< Too few pixels were assigned to this source
        MANY_PIXELS,      ///< Too many pixels were assigned to this source
        TOO_BIG,          ///< Contains pixels too far from the source center
        OVERLAP,          ///< This source overlaps with another.
        NON_POINT_SOURCE, ///< This does not appear to be a point source.
        BAD_BACKGROUND,   ///< Failed to determine a reliable background.

        ///There is a sufficient number of "better" sources than this one.
        PAST_MAX_SOURCES,

        ///There were too many saturated pixels in the source.
        MANY_SATURATED,

        NUM_DROP_REASONS, ///< How many reasons are there to drop sources.

        NOT_DROPPED = NUM_DROP_REASONS ///< The source was not dropped.
    };

    ///Human readable output of the reasons to drop sources.
    LIB_LOCAL std::ostream &operator<<(std::ostream &os,
                                       const SourceDropReason &reason);

    typedef Core::SubPixelMap GSLSubPixType;

    ///\brief Checks of whether a source should be used for PSF fitting.
    ///
    ///See get_fit_sources for a description of the arguments.
    template<class SOURCE_TYPE>
        LIB_LOCAL void check_fit_source(
            SOURCE_TYPE&              last_source,
            const Background::Source& srcbg,
            double                    max_saturated_fraction,
            unsigned                  min_pixels_per_source,
            unsigned                  max_pixels_per_source,
            double                    max_circular_aperture,
            unsigned                  min_bg_pixels,
            bool                      ignore_overlaps
        )
        {
            if(
                std::isnan(srcbg.value())
                ||
                std::isnan(srcbg.error())
                ||
                srcbg.pixels() < min_bg_pixels
            ) 
                last_source.drop(BAD_BACKGROUND);
            else if(last_source.pixel_count() > max_pixels_per_source) {
                last_source.drop(MANY_PIXELS);
            } else if(
                max_circular_aperture
                &&
                last_source.aperture() > max_circular_aperture
            ) {
                last_source.drop(TOO_BIG);
            } else if(
                last_source.saturated_pixel_count()
                >
                max_saturated_fraction * last_source.pixel_count()
            )
                last_source.drop(MANY_SATURATED);
            else if(last_source.pixel_count() < min_pixels_per_source)
                last_source.drop(FEW_PIXELS);
            else if(not ignore_overlaps && last_source.overlaps().size()) {
                for(
                    typename std::set< SOURCE_TYPE* >::const_iterator
                        overlap_iter = last_source.overlaps().begin();
                    overlap_iter != last_source.overlaps().end();
                    ++overlap_iter
                )
                    (*overlap_iter)->drop(OVERLAP);
                last_source.drop(OVERLAP);
            }
        }

    ///\brief Actually discard sources flagged as unsuitable during pixel
    ///selection.
    template<class FIT_SOURCE_TYPE>
        LIB_LOCAL void drop_unsuitable_fit_sources(
            ///The list of sources to check.
            std::list< FIT_SOURCE_TYPE * > &psf_fit_sources,

            ///The first newly extracted source in the above list.
            typename std::list< FIT_SOURCE_TYPE * >::iterator check_start,

            ///See get_fit_sources.
            std::list< FIT_SOURCE_TYPE * > &dropped_sources,

            ///An output array filled with how many sources were dropped for
            ///each drop reason.
            unsigned drop_statistics[]
        )
        {
            for(unsigned i = 0; i < NUM_DROP_REASONS; ++i)
                drop_statistics[i] = 0;

            while(check_start != psf_fit_sources.end()) {
                if((*check_start)->drop_reason() == NOT_DROPPED) {
                    ++check_start;
                } else {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Dropping source ("
                              << (*check_start)
                              << "): "
                              << (*check_start)->drop_reason()
                              << std::endl;
#endif
                    typename std::list< FIT_SOURCE_TYPE * >::iterator
                        drop_iter = check_start++;
                    ++drop_statistics[(*drop_iter)->drop_reason()];
                    (*drop_iter)->exclude_from_shape_fit();
                    if(
                        (*drop_iter)->drop_reason() == MANY_PIXELS
                        ||
                        (*drop_iter)->drop_reason() == TOO_BIG
                        ||
                        (*drop_iter)->drop_reason() == NON_POINT_SOURCE
                        ||
                        (*drop_iter)->drop_reason() == BAD_BACKGROUND
                    )
                        (*drop_iter)->exclude_from_flux_fit();
                    dropped_sources.splice(dropped_sources.end(),
                                           psf_fit_sources,
                                           drop_iter);
                }
            }
        }

    ///\brief Drop excess sources from PSF shape fitting.
    ///
    ///See get_fit_sources for a description of the arguments.
    template<class FIT_SOURCE_TYPE>
        LIB_LOCAL void trim_fit_sources(
            std::list< FIT_SOURCE_TYPE * >     &psf_fit_sources,
            unsigned                            max_sources,
            std::list< FIT_SOURCE_TYPE * >     &dropped_sources
        )
        {
            typedef typename std::list< FIT_SOURCE_TYPE * >::iterator SourceIter;
            if(psf_fit_sources.size() > max_sources) {
#ifdef TRACK_PROGRESS
                std::cerr << "Trimming the source list to " << max_sources
                          << " sources" << std::endl;
#endif
                psf_fit_sources.sort();
#ifdef TRACK_PROGRESS
                std::cerr << "Sorted by signal to noise" << std::endl;
#endif
                SourceIter first_kept = psf_fit_sources.begin();
                std::advance(first_kept, psf_fit_sources.size() - max_sources);
#ifdef TRACK_PROGRESS
                std::cerr << "Marked discarded in source assignment image"
                          << std::endl;
#endif

                for(
                    SourceIter trimmed = psf_fit_sources.begin();
                    trimmed != first_kept;
                    ++trimmed
                )
                    (*trimmed)->exclude_from_shape_fit();
                dropped_sources.splice(dropped_sources.end(),
                                       psf_fit_sources,
                                       psf_fit_sources.begin(),
                                       first_kept);
            }
        }

    ///\brief Add a newly constructed PiecewiseBidcubic source to a list.
    ///
    ///See get_fit_sources for descrption of undocumented arguments.
    LIB_LOCAL void add_new_source(
            Image<LinearSource>                     &image,
            const Core::SubPixelMap                 *subpix_map,
            const PSF::PiecewiseBicubic             &psf,
            double                                   alpha,
            double                                   max_circular_aperture,
            const std::string                       &output_fname,
            bool                                     cover_psf,

            ///The location of the source to add.
            const Core::SourceLocation &location,

            ///The backgruond to assume under the source.
            const Background::Source &srcbg,

            ///The source assignment ID of the new source.
            size_t source_assignment_id,

            ///The list to which to add the new source.
            LinearSourceList &destination
    );

    ///Select the sources to use for PSF fitting.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        LIB_PUBLIC void get_fit_sources(
            ///The image being processed. Should be a reference to the exact
            ///same variable for all sources in a single image!
            Image<FIT_SOURCE_TYPE>                        &image,

            ///The sub-pixel sensitivity map to assume. Must not be destroyed 
            ///while this object is in use.
            const Core::SubPixelMap                *subpix_map,

            ///The cental PSF coordinates of the sources in the image.
            const std::list<Core::SourceLocation>  &source_locations,

            ///The PSF to assume for the sources. Obviously parameters cannot
            ///be correctly set-up since those are being fitted, but should
            ///have the correct structure (i.e. grid for piecewis PSFs).
            const PSF_TYPE                         &psf,

            ///The minimum S/N threshold to considering a pixel above the 
            ///background
            double                                  alpha,

            ///The maximum fraction of saturated pixels for a source to be 
            ///used.
            double                                  max_saturated_fraction,

            ///The minimum number of pixels a source must have to be used.
            unsigned                                min_pixels_per_source,

            ///The maximum number of pixels allowed before excluding a
            ///source.
            unsigned                                max_pixels_per_source,

            ///The background estimate of the sources.
            Background::Measure                     &bg,

            ///The minimum number of pixels required in the background 
            ///determination.
            unsigned                                min_bg_pixels,

            ///The largest number of sources allowed in the final list.
            unsigned                                max_sources,

            ///If source pixels outside this radius are found, the source is 
            ///excluded
            double                                  max_circular_aperture,

            ///The name of the file where this source should be saved after 
            ///the fit.
            const std::string                      &output_fname,

            ///The output list of sources selected for PSF shape fitting.
            std::list< FIT_SOURCE_TYPE * >         &psf_fit_sources,

            ///The output list of sources rejected from the shape fit.
            std::list< FIT_SOURCE_TYPE * >         &dropped_sources,

            ///If true, any pixel which even partially overlaps with the PSF
            ///gets included. Otherwise, pixels are assigned by signal to
            ///noise (optionally filling up a circular aperture). This must
            ///be false for 
            bool                                    cover_psf = false,

            ///Do not drop any sources from PSF fitting (only used for zero 
            ///PSF fit at the moment).
            bool                                    do_not_drop = false,

            ///If false, any source which even partially overlaps with
            ///another is dropped from the fit.
            bool                                    ignore_overlaps = true
        )
    {
        bg.jump_to_first_source();
        std::list<Core::SourceLocation>::const_iterator
            location = source_locations.begin();
        typedef typename std::list< FIT_SOURCE_TYPE *>::iterator SourceIter;
        SourceIter first_new_source;
        
        for(
            size_t source_assignment_id = 1;
            location != source_locations.end();
            source_assignment_id++
        ) {
            Background::Source srcbg = bg();

            add_new_source(image,
                           subpix_map,
                           psf,
                           alpha,
                           max_circular_aperture,
                           output_fname,
                           cover_psf,
                           *location,
                           srcbg,
                           source_assignment_id,
                           psf_fit_sources);

            FIT_SOURCE_TYPE &last_source = *(psf_fit_sources.back());
            if ( source_assignment_id == 1 )
                first_new_source = --psf_fit_sources.end();
#ifdef TRACK_PROGRESS
            std::cerr << "Added source #"
                      << psf_fit_sources.size() 
                      << "("
                      << &last_source
                      << "), contaning "
                      << last_source.pixel_count()
                      << " pixels, with background = "
                      << last_source.background_electrons()
                      << " ("
                      << srcbg.value()
                      << ") based on "
                      << last_source.background_pixels()
                      << " pixels (" << srcbg.pixels() << ")"
                      << std::endl;
#endif
            check_fit_source(last_source,
                             srcbg,
                             max_saturated_fraction,
                             min_pixels_per_source,
                             max_pixels_per_source,
                             (cover_psf ? 0.0 : max_circular_aperture),
                             min_bg_pixels,
                             ignore_overlaps);
            ++location;
            if(location!=source_locations.end())
                if(!bg.next_source())
                    throw Error::Runtime("Smaller number of background "
                                         "measurements than sources in "
                                         "get_fit_sources!");
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Done extracting source pixels, starting source selection"
                  << std::endl;
#endif
        for(
            SourceIter src_i = first_new_source;
            src_i != psf_fit_sources.end();
            ++src_i
        )
            (*src_i)->finalize_pixels();


        unsigned drop_statistics[NUM_DROP_REASONS];
        if ( ! do_not_drop ) {
            drop_unsuitable_fit_sources(psf_fit_sources,
                                        first_new_source,
                                        dropped_sources,
                                        drop_statistics);

#ifdef VERBOSE_DEBUG
            if ( ! dropped_sources.empty() ) {
                std::cerr << "Dropped source reasons:" << std::endl;
                for(unsigned i = 0; i < NUM_DROP_REASONS; ++i)
                    std::cerr << static_cast<SourceDropReason>(i)
                              << ": " << drop_statistics[i] << std::endl;
            }
#endif

            trim_fit_sources(psf_fit_sources,
                             max_sources,
                             dropped_sources);
        }
    }

} //End FitPSF namespace.

#endif
