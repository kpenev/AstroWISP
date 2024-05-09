/**\file
 *
 * \brief Defines some of the functions needed exclusively by the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#ifndef __PSF_FITTING_H
#define __PSF_FITTING_H

#include "Config.h"
#include "IOSources.h"
#include "../Core/SharedLibraryExportMacros.h"
#include "../Background/Source.h"
#include "../Background/MeasureAnnulus.h"
#include "../Background/Zero.h"
#include "../Core/SourceLocation.h"
#include "../Core/Image.h"
#include "../Core/SubPixelMap.h"
#include "../PSF/PiecewiseBicubic.h"
#include "../IO/H5IODataTree.h"
#include "Eigen/Dense"
#include <set>
#include <iostream>

namespace FitPSF {

    template <class FIT_SOURCE_TYPE> class Image;
    template <class PSF_TYPE> class Source;

    class LinearSource;

    ///Convenience alias.
    typedef std::list<LinearSource *> LinearSourceList;

    ///Value of 1 for which we can treate a reference/pointer.
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
    LIB_LOCAL std::ostream &operator<<(
        ///The stream to write to.
        std::ostream &os,

        ///The reson for dropping the source to describe.
        const SourceDropReason &reason
    );

    ///Convenience alias.
    typedef Core::SubPixelMap GSLSubPixType;

    ///\brief Checks of whether a source should be used for PSF fitting.
    ///
    ///See select_fit_sources for a description of the arguments.
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
            else if(! ignore_overlaps && last_source.overlaps().size()) {
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

            ///See select_fit_sources.
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

    template<class FIT_SOURCE_TYPE>
        LIB_LOCAL bool order_src_pointer(
            const FIT_SOURCE_TYPE *a,
            const FIT_SOURCE_TYPE *b
        )
        {
            return *a < *b;
        }


    ///\brief Drop excess sources from PSF shape fitting.
    ///
    ///See select_fit_sources for a description of the arguments.
    template<class FIT_SOURCE_TYPE>
        LIB_LOCAL void trim_fit_sources(
            std::list< FIT_SOURCE_TYPE * >     &psf_fit_sources,
            unsigned                            max_sources,
            std::list< FIT_SOURCE_TYPE * >     &dropped_sources
        )
        {
            typedef typename std::list< FIT_SOURCE_TYPE * >::iterator SourceIter;
#ifndef NDEBUG
            typedef typename std::list< FIT_SOURCE_TYPE * >::const_iterator ConstSourceIter;
#endif

            if(psf_fit_sources.size() > max_sources) {
#ifdef TRACK_PROGRESS
                std::cerr << "Trimming the source list to " << max_sources
                          << " sources" << std::endl;
#endif
                psf_fit_sources.sort(order_src_pointer<FIT_SOURCE_TYPE>);

#ifndef NDEBUG
                ConstSourceIter previous = psf_fit_sources.begin();
                std::cerr << "Sorted sources (x, y: S/N): " << std::endl;
                for(
                        ConstSourceIter si = psf_fit_sources.begin();
                        si != psf_fit_sources.end();
                        ++si
                ) {
                    std::cerr << (*si)->x() << ", " << (*si)->y() << ": "
                        << (*si)->signal_to_noise() << " "
                        << ((**si) < (**previous))
                        << std::endl;
                    previous = si;
                }
#endif

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
    ///See select_fit_sources for descrption of undocumented arguments.
    LIB_LOCAL void add_new_source(
            ///See same name argument to select_fit_sources().
            Image<LinearSource> &image,

            ///See same name argument to select_fit_sources().
            const Core::SubPixelMap *subpix_map,

            ///See same name argument to select_fit_sources().
            const PSF::PiecewiseBicubic &psf,

            ///See same name argument to select_fit_sources().
            double alpha,

            ///See same name argument to select_fit_sources().
            double max_circular_aperture,

            ///See same name argument to select_fit_sources().
            const std::string &output_fname,

            ///See same name argument to select_fit_sources().
            bool cover_psf,

            ///The location of the source to add.
            const Core::SourceLocation &location,

            ///The values of the PSF expansion terms for this source.
            const Eigen::VectorXd &psf_terms,

            ///The backgruond to assume under the source.
            const Background::Source &srcbg,

            ///The source assignment ID of the new source.
            size_t source_assignment_id,

            ///The list to which to add the new source.
            LinearSourceList &destination
    );

    ///Select the sources to use for PSF fitting.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        LIB_PUBLIC void select_fit_sources(
            ///The image being processed. Should be a reference to the exact
            ///same variable for all sources in a single image!
            Image<FIT_SOURCE_TYPE>                 &image,

            ///The sub-pixel sensitivity map to assume. Must not be destroyed
            ///while this object is in use.
            const Core::SubPixelMap                *subpix_map,

            ///The cental PSF coordinates of the sources in the image.
            const IOSources                        &input_source_list,

            ///The PSF to assume for the sources. Obviously parameters cannot
            ///be correctly set-up since those are being fitted, but should
            ///have the correct structure (i.e. grid for piecewis PSFs).
            const PSF_TYPE                         &psf,

            ///The minimum S/N threshold to considering a pixel above the
            ///background
            double                                 alpha,

            ///The maximum fraction of saturated pixels for a source to be
            ///used.
            double                                 max_saturated_fraction,

            ///The minimum number of pixels a source must have to be used.
            unsigned                               min_pixels_per_source,

            ///The maximum number of pixels allowed before excluding a
            ///source.
            unsigned                               max_pixels_per_source,

            ///The background estimate of the sources.
            Background::Measure                    &bg,

            ///The minimum number of pixels required in the background
            ///determination.
            unsigned                               min_bg_pixels,

            ///The largest number of sources allowed in the final list.
            unsigned                               max_sources,

            ///If source pixels outside this radius are found, the source is
            ///excluded
            double                                 max_circular_aperture,

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
        std::vector<Core::SourceLocation>::const_iterator
            location = input_source_list.locations().begin();
        typedef typename std::list< FIT_SOURCE_TYPE *>::iterator SourceIter;
        SourceIter first_new_source;

        for(
            size_t source_assignment_id = 1;
            location != input_source_list.locations().end();
            source_assignment_id++
        ) {
#ifndef NDEBUG
            std::cerr << "Adding source #" << source_assignment_id
                      << ", x = " << location->x()
                      << ", y = " << location->y()
                      << std::endl;
#endif
            Background::Source srcbg = bg();

            add_new_source(
                image,
                subpix_map,
                psf,
                alpha,
                max_circular_aperture,
                input_source_list.output_fname(),
                cover_psf,
                *location,
                input_source_list.psf_terms().col(source_assignment_id - 1),
                srcbg,
                source_assignment_id,
                psf_fit_sources
            );

            FIT_SOURCE_TYPE &last_source = *(psf_fit_sources.back());
            if ( source_assignment_id == 1 )
                first_new_source = --psf_fit_sources.end();
#ifndef NDEBUG
            std::cerr << "Added source #"
                      << psf_fit_sources.size()
                      << "("
                      << &last_source
                      << ", x=" << last_source.x()
                      << ", y=" << last_source.y()
                      << "), contaning "
                      << last_source.pixel_count()
                      << " pixels, with background = "
                      << last_source.background_electrons()
                      << " ("
                      << srcbg.value()
                      << ") based on "
                      << last_source.background_pixels()
                      << " pixels (" << srcbg.pixels() << ")"
                      << " PSF map terms: " << last_source.expansion_terms()
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
            if(location!=input_source_list.locations().end())
                if(!bg.next_source())
                    throw Error::Runtime("Smaller number of background "
                                         "measurements than sources in "
                                         "select_fit_sources!");
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

            if(max_sources)
                trim_fit_sources(psf_fit_sources,
                                 max_sources,
                                 dropped_sources);
        }
    }

    ///Find the source to use for PSF fitting for a single input image.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void get_section_fit_sources(
            ///The image where fit pixels are being tracked.
            Image<FIT_SOURCE_TYPE>          &image,

            ///The configuration with which to perform PSF fitting.
            const Config                    &options,

            ///The sources in the image, location and all terms requried for
            ///PSF fitting.
            const IOSources                 &source_list,

            ///The measured background for the image sources.
            Background::Measure             &backgrounds,

            ///The sub-pixel sensitivity map to assume.
            const Core::SubPixelMap         &subpix_map,

            ///See same name argument to select_fit_sources.
            const PSF_TYPE                  &psf,

            ///See same name argument to select_fit_sources.
            std::list<FIT_SOURCE_TYPE *>    &fit_sources,

            ///See same name argument to select_fit_sources.
            std::list<FIT_SOURCE_TYPE *>    &dropped_sources
        )
        {
            double min_ston = -Core::Inf,
                   max_sat_frac = Core::Inf,
                   max_aperture = Core::Inf;
            unsigned min_pix = 0,
                     max_pix = (image.x_resolution()
                                *
                                image.y_resolution()),
                     max_src_count = std::numeric_limits<unsigned>::max();

            min_ston = options["src.min-signal-to-noise"].as<double>();
            max_sat_frac = options["src.max-sat-frac"].as<double>();
            min_pix = options["src.min-pix"].as<unsigned>();
            max_pix = options["src.max-pix"].as<unsigned>();
            max_src_count = options["src.max-count"].as<unsigned>();
            max_aperture = options["src.max-aperture"].as<double>();
#ifdef TRACK_PROGRESS
            std::cerr << "Got useful configuration for output file: ."
                      << source_list.output_fname()
                      << std::endl;
#endif

            select_fit_sources<FIT_SOURCE_TYPE, PSF_TYPE>(
                image,
                &subpix_map,
                source_list,
                psf,
                min_ston,
                max_sat_frac,
                min_pix,
                max_pix,
                backgrounds,
                options["bg.min-pix"].as<unsigned>(),
                max_src_count,
                max_aperture,
                fit_sources,
                dropped_sources,
                options.count("src.cover-bicubic-grid"),
                options["psf.model"].as<PSF::ModelType>() == PSF::ZERO
            );

            const std::vector<bool> &enabled = source_list.enabled();
            if(enabled.size())
                for(
                    typename std::list<FIT_SOURCE_TYPE *>::iterator
                    src_i = fit_sources.begin();
                    src_i != fit_sources.end();
                ) {
#ifndef NDEBUG
                    std::cerr << "Checking if source at x = "
                              << (*src_i)->x()
                              << ", y = "
                              << (*src_i)->y()
                              << " is enabled."
                              << std::endl;
#endif
                    typename std::list<FIT_SOURCE_TYPE *>::iterator
                        drop_iter = src_i++;
                    if(!enabled[(*drop_iter)->source_assignment_id() - 1]) {
                        (*drop_iter)->exclude_from_shape_fit();
                        dropped_sources.splice(dropped_sources.end(),
                                               fit_sources,
                                               drop_iter);
                    }
                }
        }

    ///True if and only if the given source is not identified by a HAT ID.
    inline bool sourceid_not_hat(const Core::SourceLocation *source)
    {
        return !(source->id().is_hatid());
    }

    ///Creates an output data tree with information common to all PSF models.
    template<class SOURCE_LIST_TYPE>
        void fill_output_data_tree_common(
            ///The sources fitted for PSF.
            const SOURCE_LIST_TYPE &fit_result,

            ///The tree to fill.
            IO::H5IODataTree &output_data_tree,

            ///The magnitude that corresponds to a flux of 1ADU.
            double mag_1adu
        )
        {
            bool hat_ids=(find_if(fit_result.begin(),
                                  fit_result.end(),
                                  sourceid_not_hat)
                          ==
                          fit_result.end());
            hat_ids = false;

            typedef std::pair< std::string, std::vector<double>* >
                DoubleKeyValue;
            typedef std::pair< std::string, std::vector<unsigned>* >
                UnsignedKeyValue;

            std::map<std::string, std::vector<double>* >
                x,
                y,
                magnitude_array,
                magnitude_error_array,
                flux_array,
                flux_error_array,
                mask_magnitude_array,
                mask_magnitude_error_array,
                mask_flux_array,
                mask_flux_error_array,
                background,
                background_error,
                chi2,
                signal_to_noise;
            std::map<std::string, std::vector<unsigned>* >
                quality_flag,
                field,
                source,
                psffit_pixels,
                background_pixels;
            std::set<std::string> output_filenames;
            std::map< std::string, std::vector<std::string>* > source_names;

            /*source names
              =new std::vector<char*>(
              hat_ids ? 0 : fit_result.size()
              );*/

            for(
                typename SOURCE_LIST_TYPE::const_iterator
                    source_i = fit_result.begin();
                source_i != fit_result.end();
                ++source_i
            ) {
                const std::string &output_fname = (*source_i)->output_filename();
                if(output_filenames.insert(output_fname).second) {
                    x.insert(DoubleKeyValue(output_fname,
                                            new std::vector<double>));
                    y.insert(DoubleKeyValue(output_fname,
                                            new std::vector<double>));
                    magnitude_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    magnitude_error_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    flux_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    flux_error_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    mask_magnitude_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    mask_magnitude_error_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    mask_flux_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    mask_flux_error_array.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    background.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    background_error.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    chi2.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    signal_to_noise.insert(
                        DoubleKeyValue(output_fname, new std::vector<double>)
                    );
                    quality_flag.insert(
                        UnsignedKeyValue(output_fname,
                                         new std::vector<unsigned>)
                    );
                    field.insert(
                        UnsignedKeyValue(output_fname,
                                         new std::vector<unsigned>)
                    );
                    source.insert(
                        UnsignedKeyValue(output_fname,
                                         new std::vector<unsigned>)
                    );
                    psffit_pixels.insert(
                        UnsignedKeyValue(output_fname,
                                         new std::vector<unsigned>)
                    );
                    background_pixels.insert(
                        UnsignedKeyValue(output_fname,
                                         new std::vector<unsigned>)
                    );
                    source_names.insert(
                        std::pair< std::string, std::vector<std::string>* >(
                            output_fname,
                            new std::vector<std::string>
                        )
                    );

                }

                x[output_fname]->push_back((*source_i)->x());
                y[output_fname]->push_back((*source_i)->y());
                magnitude_array[output_fname]->push_back(
                    magnitude((*source_i)->flux(0).value(), mag_1adu)
                );
                magnitude_error_array[output_fname]->push_back(
                    magnitude_error((*source_i)->flux(0).value(),
                                    (*source_i)->flux(0).error())
                );
                flux_array[output_fname]->push_back(
                    (*source_i)->flux(0).value()
                );
                flux_error_array[output_fname]->push_back(
                    (*source_i)->flux(0).error()
                );

                mask_magnitude_array[output_fname]->push_back(
                    magnitude((*source_i)->mask_flux().value(), mag_1adu)
                );
                mask_magnitude_error_array[output_fname]->push_back(
                    magnitude_error((*source_i)->mask_flux().value(),
                                    (*source_i)->mask_flux().error())
                );
                mask_flux_array[output_fname]->push_back(
                    (*source_i)->mask_flux().value()
                );
                mask_flux_error_array[output_fname]->push_back(
                    (*source_i)->mask_flux().error()
                );

                background[output_fname]->push_back(
                    (*source_i)->background().value()
                );
                background_error[output_fname]->push_back(
                    (*source_i)->background().error()
                );
                chi2[output_fname]->push_back((*source_i)->reduced_chi2());
                signal_to_noise[output_fname]->push_back(
                    (*source_i)->signal_to_noise()
                );
                quality_flag[output_fname]->push_back(
                    static_cast<unsigned>((*source_i)->flux(0).flag())
                );
                if(hat_ids) {
                    field[output_fname]->push_back((*source_i)->id().field());
                    source[output_fname]->push_back((*source_i)->id().source());
                } else {
                    source_names[output_fname]->push_back(
                        (*source_i)->id().str()
                    );
                }
                psffit_pixels[output_fname]->push_back(
                    (*source_i)->pixel_count()
                );
                background_pixels[output_fname]->push_back(
                    (*source_i)->background().pixels()
                );
            }

            IO::TranslateToAny< std::vector<double> > double_trans;
            IO::TranslateToAny< std::vector<unsigned> > unsigned_trans;

            typedef IO::IOTreeBase::path_type path;
            for(
                std::set<std::string>::const_iterator
                fname_i = output_filenames.begin();
                fname_i != output_filenames.end();
                ++fname_i
            ) {
                if(hat_ids) {
                    output_data_tree.put(
                        path("projsrc|srcid|field|" + *fname_i, '|'),
                        *(field[*fname_i]),
                        unsigned_trans
                    );
                    output_data_tree.put(
                        path("projsrc|srcid|source|" + *fname_i, '|'),
                        *(source[*fname_i]),
                        unsigned_trans
                    );
                } else {
                    output_data_tree.put(
                        path("projsrc|srcid|name|" + *fname_i, '|'),
                        *(source_names[*fname_i]),
                        IO::TranslateToAny< std::vector<std::string> >()
                    );
                }
                output_data_tree.put(path("projsrc|x|" + *fname_i, '|'),
                                     *(x[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("projsrc|y|" + *fname_i, '|'),
                                     *(y[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("bg|value|" + *fname_i, '|'),
                                     *(background[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("bg|error|" + *fname_i, '|'),
                                     *(background_error[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("psffit|mag|" + *fname_i, '|'),
                                     *(magnitude_array[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("psffit|mag_err|" + *fname_i, '|'),
                                     *(magnitude_error_array[*fname_i]),
                                     double_trans);
                output_data_tree.put(path("psffit|flux|" + *fname_i, '|'),
                                     *(flux_array[*fname_i]),
                                     double_trans);
                output_data_tree.put(
                    path("psffit|flux_err|" + *fname_i, '|'),
                    *(flux_error_array[*fname_i]),
                    double_trans
                );
                output_data_tree.put(
                    path("psffit|mask_mag|" + *fname_i, '|'),
                    *(mask_magnitude_array[*fname_i]),
                    double_trans
                );
                output_data_tree.put(
                    path("psffit|mask_mag_err|" + *fname_i, '|'),
                    *(mask_magnitude_error_array[*fname_i]),
                    double_trans
                );
                output_data_tree.put(
                    path("psffit|mask_flux|" + *fname_i, '|'),
                    *(mask_flux_array[*fname_i]),
                    double_trans
                );
                output_data_tree.put(
                    path("psffit|mask_flux_err|" + *fname_i, '|'),
                    *(mask_flux_error_array[*fname_i]),
                    double_trans
                );
                output_data_tree.put(path("psffit|chi2|" + *fname_i, '|'),
                                     *(chi2[*fname_i]),
                                     double_trans);
                output_data_tree.put(
                    path("psffit|sigtonoise|" + *fname_i, '|'),
                    *(signal_to_noise[*fname_i]),
                    double_trans
                );
                output_data_tree.put(path("psffit|npix|" + *fname_i, '|'),
                                     *(psffit_pixels[*fname_i]),
                                     unsigned_trans);
                output_data_tree.put(path("bg|npix|" + *fname_i, '|'),
                                     *(background_pixels[*fname_i]),
                                     unsigned_trans);
                output_data_tree.put(path("psffit|quality|" + *fname_i, '|'),
                                     *(quality_flag[*fname_i]),
                                     unsigned_trans);
                delete x[*fname_i];
                delete y[*fname_i];
                delete magnitude_array[*fname_i];
                delete magnitude_error_array[*fname_i];
                delete flux_array[*fname_i];
                delete flux_error_array[*fname_i];
                delete mask_magnitude_array[*fname_i];
                delete mask_magnitude_error_array[*fname_i];
                delete mask_flux_array[*fname_i];
                delete mask_flux_error_array[*fname_i];
                delete background[*fname_i];
                delete background_error[*fname_i];
                delete chi2[*fname_i];
                delete signal_to_noise[*fname_i];
                delete quality_flag[*fname_i];
                delete field[*fname_i];
                delete source[*fname_i];
                delete psffit_pixels[*fname_i];
                delete background_pixels[*fname_i];
                delete source_names[*fname_i];
            }
        }

    ///True if and only if the ID of s1 is less than the ID of s2.
    template<class SOURCE_TYPE>
        bool compare_source_assignment_ids(const SOURCE_TYPE *s1,
                                           const SOURCE_TYPE *s2)
        {
            return s1->source_assignment_id() < s2->source_assignment_id();
        }



} //End FitPSF namespace.

#endif
