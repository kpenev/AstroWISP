/**\file
 *
 * \brief Defines a class for working with groups of overlapping sources.
 *
 * \ingroup FitPSF
 */

#ifndef __OVERLAP_GROUPS_H
#define __OVERLAP_GROUPS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Pixel.h"
#include <set>

namespace FitPSF {

    ///Class representing a group of overlapping sources.
    template<class SOURCE_TYPE, class PSF_TYPE>
        class LIB_LOCAL OverlapGroup {
        private:
            ///Alias for the type for a list of sources.
            typedef typename std::list< SOURCE_TYPE* > SourceList;

            ///Alias for the type for a set of sources.
            typedef std::set< SOURCE_TYPE* > SourceSet;

            ///Alias for iterator type over a mutable set of sources.
            typedef typename SourceSet::iterator SourceSetIterator;

            ///Alias for iterator type over an unmutable set of sources.
            typedef typename SourceSet::const_iterator
                ConstSourceSetIterator;

            ///Alias for the type for a set of pixels.
            typedef std::set< Pixel<SOURCE_TYPE>* > PixelSet;

            ///Alias for the type for a list of pixels.
            typedef std::list< Pixel<SOURCE_TYPE>* > PixelList;

            ///Alias for iterator type over an unmutable set of pixels.
            typedef typename PixelList::const_iterator ConstPixelIter;

            ///Alias for iterator type over n mutable set of pixels.
            typedef typename PixelList::iterator PixelIter;

            ///The filename of the image this group is part of.
            unsigned __image_id;

            ///\brief The sources that belong to this group (as found in a
            ///list of dropped sources).
            SourceSet __sources;

            ///\brief The total number of pixels used for flux fitting in the
            ///group.
            ///
            ///Only valid after prepare_fitting() is called.
            unsigned __pixel_count;

            ///\brief The number of shape fitting pixels in the group.
            unsigned __shape_fit_pixel_count;

            ///The background excess values of all pixels in the group.
            Eigen::VectorXd __background_excesses;

            ///\brief Fills __background_excesses once at the end of
            ///prepare_fitting.
            void fill_background_excesses(
                ///The pixels shared by more than one sounce (assumed to
                ///coincide with the flux fitting only pixels).
                const PixelSet &shared_pixels,

                ///The common background to assume under all sources.
                double background,

                ///The variance in background.
                double background_variance
            );

            ///\brief Fills a matrix with the integrals of the PSFs of the
            ///group sources over the group pixels.
            ///
            ///Rows are pixels, columns are sources.
            template<typename PSF_INFO_TYPE>
                void fill_estimated_excess_matrix(
                    ///Sufficient information to define the PSF.
                    const PSF_INFO_TYPE &psf_info,

                    ///The matrix to fill. Resized as required.
                    Eigen::MatrixXd &estimated_excess_matrix
                );
        public:
            ///Create a new group including the given source.
            OverlapGroup(SOURCE_TYPE *first_source) :
                __image_id(first_source->image_id())
            {
                add_source(first_source);
            }

            ///Does the group include the given source (not just its ID).
            bool contains(const Core::SourceLocation &source)
            {
                return (std::find(__sources.begin(), __sources.end(), &source)
                        !=
                        __sources.end());
            }

            ///The ID of the image this group of sources belongs to
            unsigned image_id() const {return __image_id;}

            ///\brief Add the given source to the group (and all its overlaps
            ///to the IDs)
            void add_source(
                ///The source to add
                SOURCE_TYPE *source
            )
            {
                assert(source->image_id() == __image_id);

                if(__sources.insert(source).second)
                    for(
                        typename SourceSet::const_iterator
                            si = source->overlaps().begin();
                        si != source->overlaps().end();
                        ++si
                    )
                        add_source(*si);
            }

            ///The sources currently in the group (modifiable).
            SourceSet &sources() {return __sources;}

            ///The sources currently in the group (unmodifiable).
            const SourceSet &sources() const {return __sources;}

            ///\brief The number of pixels that each source in the group has.
            ///
            ///Only valid after prepare_fitting() is called.
            unsigned pixel_count() const {return __pixel_count;}

            ///Merge with the given group (and empty it).
            void merge(OverlapGroup<SOURCE_TYPE, PSF_TYPE> &group)
            {
                __sources.splice(__sources.end(), group.__sources);
            }

#ifndef NDEBUG
            ///\brief Are the overlaps of all groups sources also in the
            ///group?
            void assert_self_consistent()
            {
                for(
                    SourceSetIterator src_i = __sources.begin();
                    src_i != __sources.end();
                    ++src_i
                ) {
                    const SourceSet &overlaps = (*src_i)->overlaps();
                    for(
                        ConstSourceSetIterator overlap_i = overlaps.begin();
                        overlap_i != overlaps.end();
                        ++overlap_i
                    )
                        assert(contains(**overlap_i));
                }
            }
#endif

            ///\brief Fit for the fluxes of all sources in the group
            ///returning the sum square of the amplitude changes.
            ///
            ///Set flux() and chi2() for all sources.
            template<typename PSF_INFO_TYPE>
                double fit_fluxes(
                    ///Sufficient information to define the PSF.
                    const PSF_INFO_TYPE &psf_info
                );

            ///\brief Get ready to fit for the fluxes of the sources in the
            ///group.
            ///
            ///All sources backgrounds are set to the same value, estimated
            ///as the average over sources weighted by
            ///(number pixels)*(variance).
            ///
            ///Flux fitting only pixels get their flux fit indices changed.
            void prepare_fitting();
        }; //End OverlapGroup class.

    template<class SOURCE_TYPE, class PSF_TYPE>
        void OverlapGroup<SOURCE_TYPE, PSF_TYPE>::fill_background_excesses(
            const PixelSet &shared_pixels,
            double background,
            double background_variance
        )
        {
            assert(!std::isnan(background));
            assert(!std::isnan(background_variance));
            __background_excesses.resize(__pixel_count);
            unsigned pix_index = 0;
            for(
                ConstSourceSetIterator src_i = __sources.begin();
                src_i != __sources.end();
                ++src_i
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Group source ("
                          << *src_i
                          << ") at ("
                          << (*src_i)->x()
                          << ", "
                          << (*src_i)->y()
                          << ")"
                          << std::endl;
#endif
#ifndef NDEBUG
                unsigned start_pix_index = pix_index;
#endif
                for(
                    ConstPixelIter
                        pix_i = (*src_i)->shape_fit_pixels_begin();
                    pix_i != (*src_i)->shape_fit_pixels_end();
                    ++pix_i
                ) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Group pixel ("
                              << *pix_i
                              << ") at ("
                              << (*pix_i)->x()
                              << ", "
                              << (*pix_i)->y()
                              << ") is "
                              << ((*pix_i)->shared() ? "" : "not ")
                              << "shared and will "
                              << ((*pix_i)->shape_fit() ? "" : "not ")
                              << "participate in shape fitting with index "
                              << (*pix_i)->flux_fit_index()
                              << std::endl;
#endif
                    assert((*pix_i)->shape_fit());
                    __background_excesses[pix_index] = background_excess(
                        **pix_i,
                        background,
                        background_variance
                    );
#ifndef NDEBUG
                    if(std::isnan(__background_excesses[pix_index])) {
                        std::cerr << "Group pixel ("
                                  << *pix_i
                                  << ") at ("
                                  << (*pix_i)->x()
                                  << ", "
                                  << (*pix_i)->y()
                                  << ") is "
                                  << ((*pix_i)->shared() ? "" : "not ")
                                  << "shared and will "
                                  << ((*pix_i)->shape_fit() ? "" : "not ")
                                  << "participate in shape fitting with index "
                                  << (*pix_i)->flux_fit_index()
                                  << " value = " << (*pix_i)->measured()
                                  << " var = " << (*pix_i)->variance()
                                  << std::endl;
                    }
                    assert(!std::isnan(__background_excesses[pix_index]));
#endif
                    ++pix_index;
                }
#ifndef NDEBUG
                assert(pix_index - start_pix_index
                       ==
                       (*src_i)->shape_fit_pixel_count());
#endif
            }
            for(
                typename PixelSet::const_iterator
                    pix_i = shared_pixels.begin();
                pix_i != shared_pixels.end();
                ++pix_i
            ) {
                assert(!(*pix_i)->shape_fit());
                __background_excesses[
                    pix_index + (*pix_i)->flux_fit_index()
                ] = background_excess(
                   **pix_i,
                   background,
                   background_variance
                );
                assert(
                    !std::isnan(
                        __background_excesses[
                            pix_index + (*pix_i)->flux_fit_index()
                        ]
                    )
                );
            }
        }

    ///Output a human readable description of a group of overlapping sources.
    template<class SOURCE_TYPE, class PSF_TYPE>
        std::ostream &operator<<(
            ///The stream to output to.
            std::ostream &os,

            ///The group to describe
            const OverlapGroup<SOURCE_TYPE, PSF_TYPE> &group
        )
        {
            typedef typename std::set< SOURCE_TYPE* >::const_iterator
                SourceSetIterator;
            os << "Group: " << group.sources().size() << " sources (";
            for(SourceSetIterator s = group.sources().begin();
                s != group.sources().end();
                ++s
            ) {
                if(s != group.sources().begin()) os << ", ";
                os << (*s)->source_assignment_id();
            }
            os << ") and " << group.source_ids().size() << " source IDs (";
            for(
                std::set<unsigned>::const_iterator
                    s = group.source_ids().begin();
                s != group.source_ids().end();
                ++s
            ) {
                if(s != group.source_ids().begin()) os << ", ";
                os << *s;
            }
            os << ")";
            return os;
        }

    template<class SOURCE_TYPE, class PSF_TYPE>
        template<typename PSF_INFO_TYPE>
        void OverlapGroup<
            SOURCE_TYPE,
            PSF_TYPE
        >::fill_estimated_excess_matrix(
            const PSF_INFO_TYPE &psf_info,
            Eigen::MatrixXd &estimated_excess_matrix
        )
        {
            estimated_excess_matrix = Eigen::MatrixXd::Zero(
                __pixel_count,
                __sources.size()
            );
#ifdef VERBOSE_DEBUG
            std::cerr << "Filling excess matrix for "
                      << __sources.size()
                      << " sources, covering "
                      << __pixel_count
                      << " pixels."
                      << std::endl;
#endif
            unsigned flux_fit_pixel_count = (__pixel_count
                                             -
                                             __shape_fit_pixel_count),
                     source_column = 0,
                     shape_fit_row = 0;
            for(
                SourceSetIterator source_i = __sources.begin();
                source_i != __sources.end();
                ++source_i
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Group source ("
                          << (*source_i)
                          << ") at ("
                          << (*source_i)->x()
                          << ", "
                          << (*source_i)->y()
                          << ") has "
                          << (*source_i)->shape_fit_pixel_count()
                          << " shape fit pixels, starting from index "
                          << shape_fit_row
                          << std::endl;
#endif
                (*source_i)->fill_fluxfit_column(
                    psf_info,
                    estimated_excess_matrix.block(
                        shape_fit_row,
                        source_column,
                        (*source_i)->shape_fit_pixel_count(),
                        1
                    ),
                    estimated_excess_matrix.block(
                        __shape_fit_pixel_count,
                        source_column,
                        flux_fit_pixel_count,
                        1
                    ),
                    true,
                    false
                );

                ++source_column;
                shape_fit_row += (*source_i)->shape_fit_pixel_count();
            }
        }

    template<class SOURCE_TYPE, class PSF_TYPE>
        void OverlapGroup<SOURCE_TYPE, PSF_TYPE>::prepare_fitting()
        {
            double background = 0,
                   background_variance = 0,
                   background_norm = 0;
            PixelSet shared_pixels;
            __shape_fit_pixel_count = 0;
            bool zero_variance = false;
            for(
                SourceSetIterator source_i = __sources.begin();
                source_i != __sources.end();
                ++source_i
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Preparing group source ("
                          << *source_i
                          << ") at ("
                          << (*source_i)->x()
                          << ", "
                          << (*source_i)->y()
                          << ") for fitting!"
                          << std::endl;
#endif
                if(
                    std::isnan((*source_i)->background_electrons())
                    ||
                    std::isnan((*source_i)->background_electrons_variance())
                ) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Discarding group source ("
                              << *source_i
                              << ") at ("
                              << (*source_i)->x()
                              << ", "
                              << (*source_i)->y()
                              << ") having "
                              << (*source_i)->shape_fit_pixel_count()
                              << " shape fit and "
                              << (*source_i)->flux_fit_pixel_count()
                              << "flux fit pixels, due to bad background!"
                              << std::endl;
#endif
                    assert((*source_i)->shape_fit_pixel_count() == 0);
                    continue;
                }
                double bg_weight = 1.0;

                if((*source_i)->background_electrons_variance()) {
                    if(zero_variance)
                        bg_weight = 0.0;
                    else
                        bg_weight /=
                            (*source_i)->background_electrons_variance();
                } else zero_variance = true;

                background += (*source_i)->background_electrons() * bg_weight;
                assert(!std::isnan(background));

#ifdef VERBOSE_DEBUG
                if(background_norm > 1e10) {
                    std::cerr << "Group source ("
                              << *source_i
                              << ") at ("
                              << (*source_i)->x()
                              << ", "
                              << (*source_i)->y()
                              << ") has BG based on "
                              << (*source_i)->background_pixels()
                              << " pixels, with variance = "
                              << (*source_i)->background_electrons_variance()
                              << std::endl;
                }
#endif

                background_variance += (
                    (*source_i)->background_electrons_variance()
                    *
                    std::pow(bg_weight, 2)
                );
                assert(!std::isnan(background_variance));

                background_norm += bg_weight;
                __shape_fit_pixel_count +=
                    (*source_i)->shape_fit_pixel_count();
                for(
                    ConstPixelIter
                        pix_i = (*source_i)->flux_fit_pixels_begin();
                    pix_i != (*source_i)->flux_fit_pixels_end();
                    ++pix_i
                )
                    shared_pixels.insert(*pix_i);
            }
            assert(background_norm != 0);

            __pixel_count = __shape_fit_pixel_count + shared_pixels.size();

            background /= background_norm;
            assert(!std::isnan(background));
            assert(!std::isnan(background_variance));
            background_variance /= std::pow(background_norm, 2);
#ifdef VERBOSE_DEBUG
            if(std::isnan(background_variance)) {
                std::cerr << "NaN background variaence, from norm = "
                          << background_norm
                          << std::endl;
            }
#endif
            assert(!std::isnan(background_variance));

            for(
                SourceSetIterator source_i = __sources.begin();
                source_i != __sources.end();
                ++source_i
            ) {
                (*source_i)->set_background_electrons(background);
                (*source_i)->set_background_electrons_variance(
                    background_variance
                );
            }

            unsigned flux_fit_index = 0;
            for(
                typename PixelSet::const_iterator
                    pix_i = shared_pixels.begin();
                pix_i != shared_pixels.end();
                ++pix_i
            ) {
#ifdef VERBOSE_DEBUG
                std::cerr << "Pixel("
                          << (*pix_i)->x()
                          << ", "
                          << (*pix_i)->y()
                          << ") assigned to flux fit index "
                          << flux_fit_index
                          << std::endl;
#endif
                (*pix_i)->set_flux_fit_index(flux_fit_index++);
            }

            fill_background_excesses(shared_pixels,
                                     background,
                                     background_variance);
        }

    template<class SOURCE_TYPE, class PSF_TYPE>
        template<typename PSF_INFO_TYPE>
        double OverlapGroup<SOURCE_TYPE, PSF_TYPE>::fit_fluxes(
            const PSF_INFO_TYPE &psf_info
        )
        {
            prepare_fitting();
            typedef Eigen::JacobiSVD<
                        Eigen::MatrixXd,
                        Eigen::FullPivHouseholderQRPreconditioner
                    > SVDType;

            size_t source_count = __sources.size();
            Eigen::VectorXd amplitudes;
            double chi2;
            SVDType svd(__pixel_count,
                        source_count,
                        Eigen::ComputeFullU | Eigen::ComputeFullV);

            if(__pixel_count >= source_count) {
                Eigen::MatrixXd estimated_excess_matrix;
                fill_estimated_excess_matrix(psf_info,
                                             estimated_excess_matrix);
                assert(!std::isnan(estimated_excess_matrix.sum()));

#ifdef VERBOSE_DEBUG
                std::cerr << "Estimated excess matrix:"
                          << std::endl
                          << estimated_excess_matrix
                          << std::endl
                          << "BG excesses (size = "
                          << __background_excesses.size()
                          << "):"
                          << std::endl
                          << __background_excesses
                          << std::endl;
#endif
                svd.compute(estimated_excess_matrix);
                assert(!std::isnan(__background_excesses.sum()));
                amplitudes = svd.solve(__background_excesses);
                assert(!std::isnan(amplitudes.sum()));
                chi2 = (__background_excesses
                        -
                        estimated_excess_matrix * amplitudes).squaredNorm();
            } else {
                amplitudes = Eigen::VectorXd::Constant(source_count, Core::NaN);
                chi2 = Core::NaN;
            }
            SourceSetIterator source_iter = __sources.begin();
            double result = 0;
            for(unsigned i = 0; i < source_count; ++i) {
#ifdef DEBUG
                assert(source_iter != __sources.end());
#endif
                result += std::pow(
                    (*source_iter)->flux(0).value() - amplitudes[i],
                    2
                );
                (*source_iter)->flux(0).value() = amplitudes[i];
                (*source_iter)->flux(0).error() = (
                    __pixel_count >= source_count
                    ? (std::sqrt(chi2)
                       /
                       (
                           svd.singularValues().cwiseProduct(
                               svd.matrixV().row(i).transpose()
                           )
                       ).norm()
                       /
                       (__pixel_count - source_count))
                    : Core::NaN
                );
                (*source_iter)->chi2() = chi2;
                (*source_iter)->set_psf_amplitude(amplitudes[i]);
                (*source_iter)->set_sources_in_group(source_count);
#ifdef VERBOSE_DEBUG
                std::cerr << "Group source ("
                          << *source_iter
                          << ") at ("
                          << (*source_iter)->x()
                          << ", "
                          << (*source_iter)->y()
                          << ") flux = "
                          << (*source_iter)->flux(0).value()
                          << std::endl;
#endif
                ++source_iter;
            }
            return result;
        }

    ///Finds groups of overlapping sources among the given source list.
    template<class SOURCE_TYPE, class PSF_TYPE>
        void find_overlap_groups(
            ///The current list of shape fit sources.
            std::list< SOURCE_TYPE * > &fit_sources,

            ///The current list of sources dropped from shape fitting.
            std::list< SOURCE_TYPE * > &dropped_sources,

            ///On output this updated with new overlapping groups found. Each
            ///group is a list of iterators pointing to the  sources in the
            ///group.
            std::list< OverlapGroup<SOURCE_TYPE, PSF_TYPE> > &overlap_groups
        )
        {
            typedef std::list< OverlapGroup<SOURCE_TYPE, PSF_TYPE> >
                GroupList;
            typedef typename GroupList::iterator GroupIterator;
            typedef typename std::list< SOURCE_TYPE * > SourceList;
            typedef typename SourceList::iterator SourceListIterator;

            for(
                SourceListIterator source_i = fit_sources.begin();
                source_i != dropped_sources.end();
                ++source_i
            ) {
                if(source_i == fit_sources.end()) {
                    if(dropped_sources.size() == 0) break;
                    else source_i = dropped_sources.begin();
                }

#ifdef VERBOSE_DEBUG
                std::cerr << "Source at ("
                          << (*source_i)->x()
                          << ", "
                          << (*source_i)->y()
                          << ") overlaps with "
                          << (*source_i)->overlaps().size()
                          << " other sources."
                          << std::endl;
#endif

                if((*source_i)->overlaps().size()) {
                    unsigned image_id = (*source_i)->image_id();
                    bool new_group = true;
                    for(
                            GroupIterator
                                group_iter = overlap_groups.begin();
                            group_iter != overlap_groups.end();
                            ++group_iter
                    ) {
                        if(
                            image_id
                            ==
                            group_iter->image_id()
                        )
                            new_group = !(group_iter->contains(**source_i));
                    }
                    if(new_group)
                        overlap_groups.push_back(
                            OverlapGroup<SOURCE_TYPE, PSF_TYPE>(*source_i)
                        );
                }
            }

            for(
                SourceListIterator source_i = fit_sources.begin();
                source_i != dropped_sources.end();
                ++source_i
            ) {
                if(source_i == fit_sources.end()) {
                    if(dropped_sources.size() == 0) break;
                    else source_i = dropped_sources.begin();
                }
                (*source_i)->finalize_pixels();
            }

            for(
                GroupIterator group_iter=overlap_groups.begin();
                group_iter!=overlap_groups.end();
                ++group_iter
            ) {
#ifndef NDEBUG
                group_iter->assert_self_consistent();
#endif
                group_iter->prepare_fitting();
            }
        }

} //End FitPSF namespace.

#endif
