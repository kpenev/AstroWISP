/**\file
 *
 * \brief The implementation of the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#include "Config.h"
//#include "SDKPSFFitting.h"
#include "PiecewiseBicubic.h"
#include "IOSources.h"
#include "Source.h"
#include "../PSF/TermCalculator.h"
#include "../PSF/PiecewiseBicubic.h"
#include "../IO/SubPixHDF5File.h"
#include "../IO/FitsImage.h"
#include "../Background/Measure.h"
#include "../Background/Zero.h"
#include "../Background/MeasureAnnulus.h"
#include "../Core/PhotColumns.h"
#include "../Core/Source.h"
//#include "SubPixPhotIO.h"
#include "../Core/Error.h"
#include "../Core/SubPixelMap.h"
#include <string>
#include <list>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>

const std::string FIT_PSF_VERSION="$Revision: $";

namespace FitPSF {

    ///\brief Fit the PSF of all sources simultaneously assuming smooth 
    ///variation with position.
    ///
    ///Perform a global fit for the polynomial expansions of S, D and K over 
    ///the image using the given sources, starting from the given initial 
    ///guess if supplied.
/*    std::list<IO::OutputSDKSource> global_sdk_fit(
            ///A list of the sources in the image.
            std::list< SDKSource<Core::SubPixelMap> > &fit_sources,

            ///The fitter.
            FitPolynomialSDK<Core::SubPixelMap> &sdk_fit,

            ///The subpixel sensitivity map
            const Core::SubPixelMap &subpix_map,

            ///The fitter supports two fitting methods: Newton-Raphson and
            ///simplex. If this parameters is true the much slower simplex 
            ///method is used, otherwise Newton-Raphson.
            bool use_GSL_simplex,

            ///The gain to assume for the image.
            double gain,

            ///An initial guess for the polynomial coefficients in the 
            ///expansions of S,D and K. The order is S coefficients first, 
            ///followed by D and then K. For each variable the terms are in 
            ///increasing overall order (that is const, x, y, x^2, xy, 
            ///y^2, ...).
            const std::list<double> &initial_guess_coef=std::list<double>()
    )
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Starting global S, D, K fit using the "
                  << (use_GSL_simplex ? "simplex" : "Newton-Raphson")
                  << " method." << std::endl;
#endif
        GSLSourceIteratorType first_src = fit_sources.begin(),
        last_src = fit_sources.end();

        if(use_GSL_simplex)
            sdk_fit.gsl_fit(first_src,
                            last_src,
                            subpix_map,
                            initial_guess_coef);
        else sdk_fit.nr_fit(first_src,
                            last_src,
                            subpix_map,
                            initial_guess_coef);

        assert(sdk_fit.get_coefficients().size() % 3 == 0);
#ifdef TRACK_PROGRESS
        std::cerr << "Done with fitting." << std::endl;
#endif
        return sdk_fit.best_fit_sources(first_src,
                                        last_src,
                                        subpix_map,
                                        gain);
    }*/

    ///\brief Perform an individual fit to the PSFs of each source.
    ///
    ///The starting point can optionally be specified. The arguments are
    ///identical to those of global_sdk_fit(). See there for a description of
    ///each argument.
/*    std::list<OutputSDKSource> individual_sdk_fit(
        std::list< SDKSource<Core::SubPixelMap> > &fit_sources,
        FitPolynomialSDK<Core::SubPixelMap> &sdk_fit,
        const Core::SubPixelMap &subpix_map,
        bool use_GSL_simplex,
        double gain,
        const std::list<double> &initial_guess_coef=std::list<double>()
    )
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Starting individual S, D, K fit using the "
                  << (use_GSL_simplex ? "simplex" : "Newton-Raphson")
                  << " method." << std::endl;
        unsigned fit_progress = 0;
#endif
        std::list<OutputSDKSource> result;
        GSLSourceIteratorType src2 = fit_sources.begin(),
                              src1 = src2++;
        while(src1 != fit_sources.end()) {
#ifdef TRACK_PROGRESS
            std::cerr << "Fitting source " << fit_progress++ << "/"
                      << fit_sources.size() << std::endl;
#endif
            if(use_GSL_simplex)
                sdk_fit.gsl_fit(src1, src2, subpix_map, initial_guess_coef);
            else sdk_fit.nr_fit(src1, src2, subpix_map, initial_guess_coef);
            std::list<OutputSDKSource> fit_src=sdk_fit.best_fit_sources(
                src1++, src2++, subpix_map, gain
            );
            result.splice(result.end(), fit_src);
        }
        return result;
    }*/

    template<class SOURCE_TYPE>
        bool compare_source_assignment_ids(const SOURCE_TYPE *s1,
                                           const SOURCE_TYPE *s2)
        {
            return s1->source_assignment_id() < s2->source_assignment_id();
        }

    bool sourceid_not_hat(const Core::SourceLocation *source)
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
            double mag_1adu)
        {
            bool hat_ids=(find_if(fit_result.begin(),
                                  fit_result.end(),
                                  sourceid_not_hat)
                          ==
                          fit_result.end());

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
            std::map< std::string, std::vector< char* >* > source_names;

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
                        std::pair< std::string, std::vector< char* >* >(
                            output_fname,
                            new std::vector< char* >
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
                        new char[(*source_i)->id().str().size() + 1]
                    );
                    strcpy(source_names[output_fname]->back(),
                           (*source_i)->id().str().c_str());
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
                } else output_data_tree.put(
                    path("projsrc|srcid|name|" + *fname_i, '|'),
                    *(source_names[*fname_i]),
                    IO::TranslateToAny< std::vector<char*> >()
                );
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
            }
        }

    ///Fills an output data tree from SDK fit sources.
/*    void fill_output_data_tree_sdkfit(
        ///The sources fitted for their SDK PSF.
        const std::list<OutputSDKSource> &fit_result,

        ///The fit that was performed.
        const FitPolynomialSDK<Core::SubPixelMap> &sdk_fit,

        ///The tree to fill.
        IO::H5IODataTree &output_data_tree,

        ///Should individual S, D and K values be filled
        bool individual_sdk,

        ///The magnitude that corresponds to a flux of 1ADU/s.
        double mag_1adu)
    {
#ifdef TRACK_PROGRESS
        std::cerr << ("Filling output data tree with PSF model independent "
                      "information.") << std::endl;
#endif*/
        /*	fill_output_data_tree_common(fit_result,
            output_data_tree,
            mag_1adu);*/
/*        if(individual_sdk) {
#ifdef TRACK_PROGRESS
            std::cerr << "Adding individual S, D and K values to data tree."
                << std::endl;
#endif
            std::map<std::string, std::vector<double>* > s, d, k;
            std::set<std::string> output_filenames;

            typedef std::pair< std::string, std::vector<double>* > KeyValue;
            for(
                std::list<OutputSDKSource>::const_iterator
                    source_i = fit_result.begin();
                source_i != fit_result.end();
                ++source_i
            ) {
                const std::string
                    &output_fname = source_i->output_filename();
                if(output_filenames.insert(output_fname).second) {
                    s.insert(KeyValue(output_fname,
                                      new std::vector<double>));
                    d.insert(KeyValue(output_fname,
                                      new std::vector<double>));
                    k.insert(KeyValue(output_fname,
                                      new std::vector<double>));
                }
                s[output_fname]->push_back(source_i->psf_s());
                d[output_fname]->push_back(source_i->psf_d());
                k[output_fname]->push_back(source_i->psf_k());
            }

            typedef IO::IOTreeBase::path_type path;
            IO::TranslateToAny< std::vector<double> > double_trans;
            for(
                std::set<std::string>::const_iterator
                    fname_i = output_filenames.begin();
                fname_i != output_filenames.end();
                ++fname_i
            ) {
                output_data_tree.put(path("psffit|s|" + *fname_i, '|'),
                                     s[*fname_i],
                                     double_trans);
                output_data_tree.put(path("psffit|d|" + *fname_i, '|'),
                                     d[*fname_i],
                                     double_trans);
                output_data_tree.put(path("psffit|k|" + *fname_i, '|'),
                                     k[*fname_i],
                                     double_trans);
            }
        } else {
#ifdef TRACK_PROGRESS
            std::cerr << "Adding S, D, K PSF map to data tree." << std::endl;
#endif
            output_data_tree.put(
                "psffit.psfmap",
                sdk_fit.get_coefficients(),
                IO::TranslateToAny< std::valarray<double> >()
            );
        }
    }*/

    ///\brief For each selected and dropped source subject to PSF fitting, add
    ///the terms the PSF is allowed to depend on.
    template<class SOURCE_LIST_TYPE>
        void add_expansion_terms(
            const FitPSF::IOSources &source_list,
            const std::string &expansion_term_expression,
            SOURCE_LIST_TYPE &fit_sources,
            SOURCE_LIST_TYPE &dropped_sources
        )
    {
        std::vector<PSF::TermValarray> expansion_term_values;
        const PSF::MapVarListType &psfmap_variables = source_list.columns();
        if(expansion_term_expression != "")
            evaluate_term_expression(expansion_term_expression,
                                     psfmap_variables.begin(),
                                     psfmap_variables.end(),
                                     expansion_term_values);

        typename SOURCE_LIST_TYPE::iterator last_source = (
            dropped_sources.size() > 0 ? dropped_sources.end()
                                       : fit_sources.end()
        );
        for(
            typename SOURCE_LIST_TYPE::iterator src_i = fit_sources.begin();
            src_i != last_source;
            ++src_i
        ) {
            if(src_i == fit_sources.end()) src_i = dropped_sources.begin();

            Eigen::VectorXd &expansion_terms = (*src_i)->expansion_terms();
            expansion_terms.resize(expansion_term_expression == ""
                                   ? 1
                                   : expansion_term_values.size());
            for(
                unsigned term_i = 0;
                term_i < expansion_terms.size();
                ++term_i
            )
                if(expansion_term_expression == "") {
                    expansion_terms[term_i] = 1;
                } else {
                    expansion_terms[term_i] = expansion_term_values[term_i][
                        (*src_i)->source_assignment_id() - 1
                    ];
                }
#ifdef VERBOSE_DEBUG
            std::cerr << "Source("
                      << (*src_i)->x()
                      << ", "
                      << (*src_i)->y()
                      << ") terms("
                      << (*src_i)->expansion_terms().size()
                      << "):";
            for(
                unsigned term_i = 0;
                term_i < (*src_i)->expansion_terms().size();
                ++term_i
            )
                std::cerr << " " << (*src_i)->expansion_terms()[term_i];
            std::cerr << std::endl;
#endif
        }
    }

    ///Find the source to use for PSF fitting for a single input image.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        void get_section_fit_sources(
            FitPSF::Image<FIT_SOURCE_TYPE>      &image,
            const FitPSF::Config                &options,
            const FitPSF::IOSources             &source_list,
            const Core::SubPixelMap             &subpix_map,
            const PSF_TYPE                      &psf,
            std::list<FIT_SOURCE_TYPE *>        &fit_sources,
            std::list<FIT_SOURCE_TYPE *>        &dropped_sources)
        {
            double min_ston = -Core::Inf,
                   max_sat_frac = Core::Inf,
                   max_aperture = Core::Inf;
            unsigned min_pix = 0,
                     max_pix = (image.x_resolution()
                                *
                                image.y_resolution()),
                     max_src_count = std::numeric_limits<unsigned>::max();

            if(options["psf.model"].as<PSF::ModelType>() != PSF::ZERO) {
                min_ston = options["src.min-signal-to-noise"].as<double>();
                max_sat_frac = options["src.max-sat-frac"].as<double>();
                min_pix = options["src.min-pix"].as<unsigned>();
                max_pix = options["src.max-pix"].as<unsigned>();
                max_src_count = options["src.max-count"].as<unsigned>();
                max_aperture = options["src.max-aperture"].as<double>();
            }

            Background::Measure *backgrounds;
            if(options["bg.zero"].as<bool>()) {
                backgrounds = new Background::Zero(
                    source_list.locations().size()
                );
            } else {
                const Background::Annulus& background_annulus =
                    options["bg.annulus"].as<Background::Annulus>();
                backgrounds = new Background::MeasureAnnulus(
                    background_annulus.inner_radius(),
                    background_annulus.outer_radius(),
                    background_annulus.inner_radius(),
                    image,
                    source_list.locations()
                );
            }

            get_fit_sources<FIT_SOURCE_TYPE, PSF_TYPE>(
                image,
                &subpix_map,
                source_list.locations(),
                psf,
                min_ston,
                max_sat_frac,
                min_pix,
                max_pix,
                *backgrounds,
                options["bg.min-pix"].as<unsigned>(),
                max_src_count,
                max_aperture,
                source_list.output_fname(),
                fit_sources,
                dropped_sources,
                options.count("src.cover-bicubic-grid"),
                options["psf.model"].as<PSF::ModelType>() == PSF::ZERO
            );

            delete backgrounds;

        }

    ///Find the sources on which PSF fitting will be performed.
    template<class FIT_SOURCE_TYPE, class PSF_TYPE>
        std::list< FitPSF::Image<FIT_SOURCE_TYPE>* > prepare_fit_sources(
            const FitPSF::Config &options,
            std::list<FIT_SOURCE_TYPE *> &fit_sources,
            std::list<FIT_SOURCE_TYPE *> &dropped_sources,
            unsigned &x_resolution,
            unsigned &y_resolution,
            const Core::SubPixelMap &subpix_map,
            const PSF_TYPE &psf,
            IO::H5IODataTree &output_data_tree
        )
        {
            std::string source_list_fname = 
                options["io.source-list"].as<std::string>();

            std::istream *source_list_stream;
            if(source_list_fname != "")
                source_list_stream = new std::ifstream(
                    source_list_fname.c_str()
                );
            else source_list_stream = &std::cin;

            std::list< FitPSF::Image<FIT_SOURCE_TYPE>* > fit_images;

            x_resolution = y_resolution = 0;
            while(true) {
                FitPSF::IOSources source_list(
                    *source_list_stream,
                    options["io.input-columns"].as<Core::StringList>()
                );

                FitPSF::Image<FIT_SOURCE_TYPE> *image;
                
                if(options["io.expect-error-hdu"].as<unsigned>() > 0) {
                    image = new FitPSF::Image<FIT_SOURCE_TYPE>(
                        source_list.fits_fname(),
                        0,
                        options["io.expect-error-hdu"].as<unsigned>()
                    );
#ifdef DEBUG
                    std::cerr << "Read image with error HDU." << std::endl;
#endif
                    assert(image->has_errors());
                } else {
                    image = new FitPSF::Image<FIT_SOURCE_TYPE>(
                        source_list.fits_fname(),
                        0,
                        options["gain"].as<double>()
                    );
#ifdef DEBUG
                    std::cerr << "Read image without error HDU." << std::endl;
#endif
                    assert(!image->has_errors());
                }

                if(x_resolution == 0 && y_resolution == 0) {
                    x_resolution = image->x_resolution();
                    y_resolution = image->y_resolution();
                } else if(
                    x_resolution != image->x_resolution()
                    ||
                    y_resolution != image->y_resolution()
                ) {
                    std::ostringstream message;
                    message << "Image " << source_list.fits_fname()
                            << " has a different resolution ("
                            << image->x_resolution()
                            << " x "
                            << image->y_resolution()
                            << ") than previous images ("
                            << x_resolution << " x " << y_resolution;
                    throw Error::FitsImage(message.str());
                }

                std::list<FIT_SOURCE_TYPE *> section_fit_sources,
                                             section_dropped_sources;
                get_section_fit_sources<FIT_SOURCE_TYPE, PSF_TYPE>(
                    *image,
                    options,
                    source_list,
                    subpix_map,
                    psf,
                    section_fit_sources,
                    section_dropped_sources
                );
                add_expansion_terms(source_list,
                                    options["psf.terms"].as<std::string>(),
                                    section_fit_sources,
                                    section_dropped_sources);
                if(options["psf.terms"].as<std::string>() != "") {
                    typedef IO::IOTreeBase::path_type path;
                    output_data_tree.put(
                        path(
                            "psffit|variables|" + source_list.output_fname(),
                            '|'
                        ),
                        source_list.columns(),
                        IO::TranslateToAny<PSF::MapVarListType>()
                    );
                }
                fit_sources.splice(fit_sources.end(), section_fit_sources);
                dropped_sources.splice(dropped_sources.end(),
                                       section_dropped_sources);

                fit_images.push_back(image);

                if(source_list.last()) break;
            }
            if(source_list_fname != "") {
                delete source_list_stream;
            }
            return fit_images;
        }

    ///\brief Does not actually fit for the PSF but sets up a constant one.
    bool zero_fit(const FitPSF::Config &options,
                  const Core::SubPixelMap &subpix_map,
                  IO::H5IODataTree &output_data_tree)
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Extracting source pixels for zero PSF fitting."
                  << std::endl;
#endif
/*        const PSF::Grid &options_grid = (
            options["psf.bicubic.grid"].as<PSF::Grid>()
        );*/
        PSF::Grid grid;
        grid.x_grid.resize(2, 0.0);
        grid.y_grid.resize(2, 0.0);
/*        grid.x_grid[0] = options_grid.x_grid.front();
        grid.x_grid[1] = options_grid.x_grid.back();
        grid.y_grid[0] = options_grid.y_grid.front();
        grid.y_grid[1] = options_grid.y_grid.back();*/

        PSF::PiecewiseBicubic psf(grid.x_grid.begin(),
                                  grid.x_grid.end(),
                                  grid.y_grid.begin(),
                                  grid.y_grid.end());
        std::vector<double> zeros(grid.x_grid.size() * grid.y_grid.size(),
                                  0);
        psf.set_values(zeros.begin(), zeros.begin(),
                       zeros.begin(), zeros.begin());

        unsigned x_resolution, y_resolution;
        LinearSourceList fit_sources, dropped_sources;
        std::list< FitPSF::Image<LinearSource>* > fit_images = 
            prepare_fit_sources<LinearSource, PSF::PiecewiseBicubic>(
                options,
                fit_sources,
                dropped_sources,
                x_resolution,
                y_resolution,
                subpix_map,
                psf,
                output_data_tree
            );

        if ( ! dropped_sources.empty() ) {
            throw Error::Runtime("Not all sources in the image were accepted "
                                 "for zero PSF fitting." );
        }
        for(
            LinearSourceList::iterator si=fit_sources.begin();
            si!=fit_sources.end();
            ++si
        ) {
            (*si)->flux(0).value() = 1.0;
            (*si)->flux(0).error() = 0.0;
            (*si)->flux(0).flag() = Core::GOOD;
            (*si)->chi2() = Core::NaN;
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Filling output data tree." << std::endl;
#endif
        fill_output_data_tree_common(fit_sources,
                                     output_data_tree,
                                     options["magnitude-1adu"].as<double>());

        output_data_tree.put("psffit.grid",
                             std::string(grid),
                             IO::translate_string);

        Eigen::VectorXd best_fit_coef;
        output_data_tree.put("psffit.psfmap",
                             best_fit_coef,
                             IO::TranslateToAny<Eigen::VectorXd>());

        for(
            std::list< FitPSF::Image<LinearSource>* >::iterator
                img_i = fit_images.begin();
            img_i != fit_images.end();
            ++img_i
        )
            delete *img_i;
        return true;
    }

    ///\brief Fit an elliptical gaussian PSF model according to the command 
    ///line options.
/*    bool sdk_fit(const FitPSF::Config &options,
                 const Core::SubPixelMap &subpix_map,
                 IO::H5IODataTree &output_data_tree)
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Extracting source pixels for SDK fitting."
                  << std::endl;
#endif
        std::list< SDKSource<Core::SubPixelMap> > fit_sources,
                                                  dropped_sources;
        unsigned x_resolution, y_resolution;
        prepare_fit_sources(options,
                            fit_sources, dropped_sources,
                            x_resolution, y_resolution,
                            output_data_tree);


        bool use_simplex=(options.count("psf.sdk.use-simplex")!=0
                          &&
                          options["psf.sdk.use-simplex"].as<bool>());
        if(!use_simplex)
            for(
                std::list< SDKSource<Core::SubPixelMap> >::iterator
                    si = fit_sources.begin();
                si != fit_sources.end();
                ++si
            ) {
                si->enable_first_deriv();
                si->enable_second_deriv();
            }
        unsigned num_terms = fit_sources.front().expansion_terms().size();
        FitPolynomialSDK<Core::SubPixelMap> sdk_fit(
            options["psf.sdk.minS"].as<double>(),
            options["psf.sdk.maxS"].as<double>(),
            options["psf.sdk.fit-tolerance"].as<double>(),
            (
                num_terms > 0
                ? options["psf.max-chi2"].as<double>()
                : Core::Inf
            )
        );
        std::list<double> initial_guess_coef =
            read_sdk_coef(options["io.initial-guess"].as<std::string>());
        std::list<OutputSDKSource> output_sources;
        if(options["psf.terms"].as<std::string>() == "")
            output_sources = individual_sdk_fit(fit_sources,
                                              sdk_fit,
                                              subpix_map,
                                              use_simplex,
                                              options["gain"].as<double>(),
                                              initial_guess_coef);
        else output_sources = global_sdk_fit(fit_sources,
                                             sdk_fit,
                                             subpix_map,
                                             use_simplex,
                                             options["gain"].as<double>(),
                                             initial_guess_coef);
#ifdef TRACK_PROGRESS
        std::cerr << "Adding fit information to output data tree." << std::endl;
#endif
        fill_output_data_tree_sdkfit(
            output_sources,
            sdk_fit,
            output_data_tree,
            options["psf.terms"].as<std::string>() == "",
            options["magnitude-1adu"].as<double>()
        );
        return true;
    }*/

    ///\brief Fit a piecewise bicubic PSF model according to the command line
    ///options.
    bool piecewise_bicubic_fit(const FitPSF::Config &options,
                               const Core::SubPixelMap &subpix_map,
                               IO::H5IODataTree &output_data_tree)
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Starting piecewise bicubic fit."
                  << std::endl;
#endif

        unsigned x_resolution,
                 y_resolution;
        LinearSourceList fit_sources,
                         dropped_sources;

        const PSF::Grid& grid = (
            options["psf.bicubic.grid"].as<PSF::Grid>()
        );
        PSF::PiecewiseBicubic psf(grid.x_grid.begin(),
                                  grid.x_grid.end(),
                                  grid.y_grid.begin(),
                                  grid.y_grid.end());
        std::vector<double> zeros(grid.x_grid.size() * grid.y_grid.size(),
                                  0);
        psf.set_values(zeros.begin(), zeros.begin(),
                       zeros.begin(), zeros.begin());

        std::list< FitPSF::Image<LinearSource>* > fit_images = 
            prepare_fit_sources<LinearSource, PSF::PiecewiseBicubic>(
                options,
                fit_sources,
                dropped_sources,
                x_resolution,
                y_resolution,
                subpix_map,
                psf,
                output_data_tree
            );

        double gain = options["gain"].as<double>();

        if ( fit_sources.empty() ) {
            throw Error::Runtime( "no valid source found in the image" );
        }

        Eigen::VectorXd best_fit_coef;
        bool converged;
#ifdef TRACK_PROGRESS
        std::cerr << "Starting fit." << std::endl;
#endif
        bool ignore_dropped = (options.count("psf.ignore-dropped") != 0
                               &&
                               options["psf.ignore-dropped"].as<bool>());
        LinearSourceList empty_source_list;
        if(options["io.initial-guess"].as<std::string>() != "") {
            IO::H5IODataTree initial_guess_data;
            IO::SubPixHDF5File initial_guess_file(
                options["io.initial-guess"].as<std::string>().c_str(),
                H5F_ACC_RDONLY
            );
            initial_guess_file.read(
                PSF::PiecewiseBicubicMap::required_data_tree_quantities().begin(),
                PSF::PiecewiseBicubicMap::required_data_tree_quantities().end(),
                initial_guess_data
            );
            PSF::PiecewiseBicubicMap psf_map(initial_guess_data);
            converged = fit_piecewise_bicubic_psf(
                fit_sources,
                (ignore_dropped ? empty_source_list : dropped_sources),
                psf_map,
                grid.x_grid,
                grid.y_grid,
                subpix_map,
                options["psf.bicubic.max-abs-amplitude-change"].as<double>(),
                options["psf.bicubic.max-rel-amplitude-change"].as<double>(),
                options["psf.max-chi2"].as<double>(),
                options["psf.bicubic.pixrej"].as<double>(),
                options["psf.min-convergence-rate"].as<double>(),
                options["psf.max-iterations"].as<int>(),
                options["psf.bicubic.smoothing"].as<double>(),
                best_fit_coef,
                gain
            );
        } else {
            converged=fit_piecewise_bicubic_psf(
                    fit_sources,
                    (ignore_dropped ? empty_source_list : dropped_sources),
                    gain,
                    grid.x_grid,
                    grid.y_grid,
                    subpix_map,
                    options[
                        "psf.bicubic.max-abs-amplitude-change"
                    ].as<double>(),
                    options[
                        "psf.bicubic.max-rel-amplitude-change"
                    ].as<double>(),
                    options["psf.max-chi2"].as<double>(),
                    options["psf.bicubic.pixrej"].as<double>(),
                    options["psf.min-convergence-rate"].as<double>(),
                    options["psf.max-iterations"].as<int>(),
                    options["psf.bicubic.smoothing"].as<double>(),
                    best_fit_coef
            );
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Adding dropped sources." << std::endl;
#endif
        fit_sources.splice(fit_sources.end(), dropped_sources);
        fit_sources.sort(compare_source_assignment_ids<LinearSource>);
        PSF::PiecewiseBicubicMap psf_map(best_fit_coef,
                                         grid.x_grid,
                                         grid.y_grid);
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            PSF::PiecewiseBicubic *psf = psf_map((*si)->expansion_terms());
            (*si)->calculate_mask_flux(*psf);
            delete psf;
            (*si)->flux(0).flag() = (*si)->quality_flag();
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Filling output data tree." << std::endl;
#endif
        fill_output_data_tree_common(fit_sources,
                                     output_data_tree,
                                     options["magnitude-1adu"].as<double>());
        output_data_tree.put("psffit.psfmap",
                             best_fit_coef,
                             IO::TranslateToAny< Eigen::VectorXd >());

        for(
            std::list< FitPSF::Image<LinearSource>* >::iterator
                img_i = fit_images.begin();
            img_i != fit_images.end();
            ++img_i
        )
            delete *img_i;

        return converged;
    }

    ///Fills the sub-pixel map to use for PSF fitting.
    void fill_subpix_map(
            ///The parsed command line options.
            const FitPSF::Config &options,

            ///The sub-pixel map to fill.
            Core::SubPixelMap &subpix_map)
    {
        if(options.executable() == "fitprf") {
            subpix_map.set_resolution(0, 0);
        } else {
            assert(options.executable() == "fitpsf");
            if(options["io.subpix"].as<std::string>() != "") {
                IO::FitsImage<double>
                    subpix_map_data(options["io.subpix"].as<std::string>());
                unsigned xres = subpix_map_data.x_resolution(),
                         yres = subpix_map_data.y_resolution();
                subpix_map.set_resolution(xres, yres);
                for(unsigned y = 0; y < yres; y++)
                    for(unsigned x = 0; x < xres; x++)
                        subpix_map(x, y) = subpix_map_data(x, y);
            } else {
                subpix_map.set_resolution(1, 1);
                subpix_map(0, 0) = 1;
            }
        }
    }

    void output(const IO::H5IODataTree &output_data_tree)
    {
#ifdef TRACK_PROGRESS
        std::cerr << "Outputting data tree:" << std::endl
                  << output_data_tree
                  << std::endl;
#endif
        const IO::IOTreeBase &psffit_var_tree =
            output_data_tree.get_child("bg.value");
        for(
            IO::IOTreeBase::const_iterator var_i = psffit_var_tree.begin();
            var_i != psffit_var_tree.end();
            ++var_i
        ) {
#ifdef DEBUG
            std::cout << "Output file name: " << var_i->first << std::endl;
#endif
            IO::SubPixHDF5File *file;
            try {
                file = new IO::SubPixHDF5File(var_i->first.c_str(),
                                              H5F_ACC_RDWR);
            } catch(H5::FileIException) {
                try {
                    file=new IO::SubPixHDF5File(var_i->first.c_str(),
                                                H5F_ACC_TRUNC);
                } catch(H5::FileIException &ex) {
                    ex.printErrorStack();
                    continue;
                }
            }
            file->write(output_data_tree, false);
            file->close();
            delete file;
        }
    }

} //End FitPSF namespace.

///Perform the PSF fitting according to the command line options.
int main(int argc, char *argv[])
{
    std::cerr.setf(std::ios::scientific);
    std::cerr.precision(16);
#ifndef DEBUG
	try {
#endif
        FitPSF::Config options(argc, argv);
		if(!options.proceed()) return 1;
#ifdef TRACK_PROGRESS
		std::cerr << "Parsed command line." << std::endl;
#endif

        Core::SubPixelMap subpix_map(0, 0, "");
        FitPSF::fill_subpix_map(options, subpix_map);

#ifdef TRACK_PROGRESS
		std::cerr << "Created sub-pixel map." << std::endl;
#endif

        IO::H5IODataTree output_data_tree(argc,
                                          argv,
                                          FIT_PSF_VERSION,
                                          options);

#ifdef TRACK_PROGRESS
		std::cerr << "Constructed output data tree." << std::endl;
#endif

		bool converged;
        if (options["psf.model"].as<PSF::ModelType>() == PSF::ZERO)
            converged = FitPSF::zero_fit(options,
                                         subpix_map,
                                         output_data_tree);
/*        else if(options["psf.model"].as<PSF::ModelType>() == PSF::SDK)
			converged = sdk_fit(options, subpix_map, output_data_tree);*/
		else 
            converged = piecewise_bicubic_fit(options,
                                              subpix_map,
                                              output_data_tree);

#ifdef TRACK_PROGRESS
		std::cerr << "Writing output HDF5 file." << std::endl;
#endif

        FitPSF::output(output_data_tree);
        if(converged) return 0;
        else
            throw Error::Fitting(
                "PSF fitting failed to converge! Accepting last iteration."
            );

        return 0;
#ifdef TRACK_PROGRESS
		std::cerr << "Done" << std::endl;
#endif

#ifndef DEBUG
	} catch(Error::General &ex) {
		std::cerr << ex.what() << ":" << ex.get_message() << std::endl;
		return 2;
    } catch(H5::Exception &ex) {
        std::cerr 
            << ex.getFuncName()
            << ": "
            << ex.getDetailMsg()
            << std::endl;
        return 3;
	} catch(std::exception &ex) {
		std::cerr << ex.what() << std::endl;
		return 4;
	}
#endif
}
