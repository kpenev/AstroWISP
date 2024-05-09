/**\file
 *
 * \brief Define the C-interface functions for the FitPSF library.
 *
 * \ingroup FitPSF
 */

#include "CInterface.h"
#include "Config.h"
#include "Image.h"
#include "LinearSource.h"
#include "PiecewiseBicubic.h"
#include "../IO/CommandLineConfig.h"
#include <cstdarg>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

FittingConfiguration *create_psffit_configuration()
{
    FittingConfiguration* result = reinterpret_cast<FittingConfiguration*>(
        new FitPSF::Config(0, NULL)
    );
#ifdef VERBOSE_DEBUG
    std::cerr << "craeted fitpsf config at: " << result << std::endl;
#endif
    return result;
}

void destroy_psffit_configuration(
    FittingConfiguration *configuration
)
{
#ifdef VERBOSE_DEBUG
    std::cerr << "destroying fitpsf config at: " << configuration << std::endl;
#endif
    delete reinterpret_cast<FitPSF::Config*>(configuration);
}

void update_psffit_configuration(bool prffit,
                                 FittingConfiguration *target_configuration,
                                 ...)
{
    FitPSF::Config *configuration =
        reinterpret_cast<FitPSF::Config*>(target_configuration);

    va_list options;
    va_start(options, target_configuration);

    IO::update_configuration(*configuration,
                             (prffit ? "fitprf" : "fitpsf"),
                             options);

    va_end(options);
}

///Create a list of all fitting sources (assign pixels, filter etc.).
void prepare_fit_sources(
    ///The configuration for PSF fitting.
    const FitPSF::Config &configuration,

    ///The PSF fitting images for simultaneous fits.
    std::vector< FitPSF::Image<FitPSF::LinearSource> > &fit_images,

    ///See same name argument to piecewise_bicubic_fit()
    char ***source_ids,

    ///See same name argument to piecewise_bicubic_fit()
    double **source_coordinates,

    ///See same name argument to piecewise_bicubic_fit()
    double **psf_terms,

    ///See same name argument to piecewise_bicubic_fit()
    bool **enabled,

    ///See same name argument to piecewise_bicubic_fit()
    unsigned long *number_sources,

    ///See same name argument to piecewise_bicubic_fit()
    unsigned long number_terms,

    ///The measured background for the sources in each image, indexed by the
    ///image index.
    BackgroundMeasureAnnulus **backgrounds,

    ///The list to add the newly created sources suitable for participating in
    ///the shape fit.
    FitPSF::LinearSourceList &fit_sources,

    ///The list to add the newly created sources if they are not suitable for
    ///participating in the shape fit.
    FitPSF::LinearSourceList &dropped_sources,

    ///The sub-pixel sensitivity map to assume.
    const Core::SubPixelMap &subpix_map,

    ///An instance of the PSF to use for PSF fitting, with properly defined
    ///grid, obviously no coefficients.
    const PSF::PiecewiseBicubic &psf,

    ///The object to add result data to (e.g. PSF map terms).
    IO::H5IODataTree &output_data_tree
)
{
    for(
        unsigned long image_index = 0;
        image_index < fit_images.size();
        ++image_index
    ) {
#ifdef TRACK_PROGRESS
        std::cerr << "Extracting sources from image "
                  << image_index
                  << std::endl;
#endif
        std::ostringstream image_index_stream;
        image_index_stream << image_index;

        FitPSF::IOSources image_sources(image_index_stream.str().c_str(),
                                        source_ids[image_index],
                                        source_coordinates[image_index],
                                        psf_terms[image_index],
                                        enabled[image_index],
                                        number_sources[image_index],
                                        number_terms);
#ifndef NDEBUG
        std::cerr << "List contains "
                  << image_sources.locations().size()
                  << " sources";

        std::cerr << " located at: " << std::endl;
        for(unsigned i = 0; i < number_sources[image_index]; ++i)
            std::cerr << std::setw(25)
                      << source_coordinates[image_index][2 * i]
                      << ", "
                      << std::setw(25)
                      << source_coordinates[image_index][2 * i + 1]
                      << std::endl;

        std::cerr << std::endl;
#endif

        FitPSF::LinearSourceList section_fit_sources,
                                 section_dropped_sources;
        FitPSF::get_section_fit_sources<FitPSF::LinearSource,
                                        PSF::PiecewiseBicubic>(
            fit_images[image_index],
            configuration,
            image_sources,
            *reinterpret_cast<Background::MeasureAnnulus*>(
                backgrounds[image_index]
            ),
            subpix_map,
            psf,
            section_fit_sources,
            section_dropped_sources
        );
#ifdef TRACK_PROGRESS
        std::cerr << "Extracted "
                  << section_fit_sources.size()
                  << " sources."
                  << std::endl;
#endif

#ifndef NDEBUG
        std::cerr << "Selected section source locations and S/N: " << std::endl;
#endif
        for(
                FitPSF::LinearSourceList::const_iterator
                si = section_fit_sources.begin();
                si != section_dropped_sources.end();
                ++si
        ) {
            if(si == section_fit_sources.end()) {
                si = section_dropped_sources.begin();
#ifndef NDEBUG
                std::cerr << "Dropped sources: " << std::endl;
#endif
            }
            if(si != section_dropped_sources.end()) {
                (*si)->finalize_pixels();
#ifndef NDEBUG
                std::cerr << (*si)->x() << ", " << (*si)->y()
                          << ": " << (*si)->signal_to_noise() << std::endl;
#endif
            }
        }


        output_data_tree.put(
            "psffit.terms." + image_index_stream.str(),
            std::vector<double>(
                image_sources.psf_terms().data(),
                (
                    image_sources.psf_terms().data()
                    +
                    number_sources[image_index] * number_terms
                )
            ),
            IO::TranslateToAny< std::vector<double> >()
        );

#ifndef NDEBUG
        std::cerr << "Trying to read back psffit terms. Node name: "
                  << "psffit.terms." + image_index_stream.str()
                  << std::endl;
        output_data_tree.get< std::vector<double> >(
            "psffit.terms." + image_index_stream.str(),
            std::vector<double>(),
            IO::TranslateToAny< std::vector<double> >()
        );
        std::cerr << "Finished reading back psffit terms." << std::endl;
#endif

#ifdef TRACK_PROGRESS
        std::cerr << "Added PSF fit terms to result tree." << std::endl;
#endif
        fit_sources.splice(fit_sources.end(), section_fit_sources);
        dropped_sources.splice(dropped_sources.end(), section_dropped_sources);
#ifdef TRACK_PROGRESS
        std::cerr << "Integrated into global source lists:"
                  << fit_sources.size() << " fit sources and "
                  << dropped_sources.size() << " dropped sources"
                  << std::endl;
#endif
    }
}

///Get the variables that participate in the PSF map from a data tree.
LIB_PUBLIC bool local_get_psf_map_variables(
    ///The data tree containing the variables to extract.
    H5IODataTree *output_data_tree,

    ///The index of the image for wich to get the variables within the list of
    ///input images.
    unsigned image_index,

    ///An array to fill with the values of the variables. All values of the
    ///first variable come first, followed by all values of the second variable
    ///etc.
    double *column_data
)
{
    IO::H5IODataTree *real_output_data_tree =
        reinterpret_cast<IO::H5IODataTree*>(output_data_tree);

    std::ostringstream tree_path;
    tree_path << "psffit.variables." << image_index;
#ifdef VERBOSE_DEBUG
    std::cerr << "Getting PSF map variables from tree at "
              << real_output_data_tree
              << " with path "
              << tree_path.str()
              << std::endl;
#endif

    const PSF::MapVarListType &variables =
        real_output_data_tree->get<PSF::MapVarListType>(
            std::string("psffit.variables.0"),
            PSF::MapVarListType(),
            IO::TranslateToAny<PSF::MapVarListType>()
        );
    if(variables.size() == 0) {
#ifdef VERBOSE_DEBUG
        std::cerr << "Empty varibales entry in result tree!" << std::endl;
#endif
        return false;
    }
#ifdef VERBOSE_DEBUG
    else
        std::cerr << "Finished reading back psffit variables with zise:"
                  << variables.size() << "x" << variables.begin()->second.size()
                  << std::endl;
#endif

    double *destination = column_data;
    for(
        PSF::MapVarListType::const_iterator var_i = variables.begin();
        var_i != variables.end();
        ++var_i
    ) {
        const double *start = &(var_i->second[0]),
                     *end = start + var_i->second.size();

        std::copy(start, end, destination);
        destination += var_i->second.size();
    }
    return true;
}

bool piecewise_bicubic_fit(double **pixel_values,
                           double **pixel_errors,
                           char **pixel_masks,
                           unsigned long number_images,
                           unsigned long image_x_resolution,
                           unsigned long image_y_resolution,
                           char ***source_ids,
                           double **source_coordinates,
                           double **psf_terms,
                           bool **enabled,
                           unsigned long *number_sources,
                           unsigned long number_terms,
                           BackgroundMeasureAnnulus** backgrounds,
                           FittingConfiguration *configuration,
                           double *subpix_sensitivities,
                           unsigned long subpix_x_resolution,
                           unsigned long subpix_y_resolution,
                           H5IODataTree *output_data_tree)
{
#ifdef TRACK_PROGRESS
    std::cerr << "Starting piecewise bicubic fit." << std::endl;
#endif

#ifndef NDEBUG
    std::cerr << "Source locations and PSF terms: " << std::endl;
    for(unsigned image_i = 0; image_i < number_images; ++image_i) {
        std::cerr << "Image " << image_i << std::endl;
        for(
                unsigned source_i = 0;
                source_i < number_sources[image_i];
                ++source_i
        ) {
            std::cerr << "\t"
                      << std::setw(25)
                      << source_coordinates[image_i][2 * source_i]
                      << ", "
                      << std::setw(25)
                      << source_coordinates[image_i][2 * source_i + 1];
            for(
                unsigned term_i = 0;
                term_i < number_terms;
                ++term_i
            )
                std::cerr
                    << ", "
                    << std::setw(25)
                    << psf_terms[image_i][number_terms * source_i + term_i];
            std::cerr << std::endl;
        }
    }
#endif
    Core::SubPixelMap subpix_map(subpix_sensitivities,
                                 subpix_x_resolution,
                                 subpix_y_resolution);

#ifdef TRACK_PROGRESS
    std::cerr << "Created subpixel map." << std::endl;
#endif
    std::vector< FitPSF::Image<FitPSF::LinearSource> >
        fit_images(number_images);
    for(
        unsigned long image_index = 0;
        image_index < number_images;
        ++image_index
    )
        fit_images[image_index].wrap(
            pixel_values[image_index],
            pixel_masks[image_index],
            image_x_resolution,
            image_y_resolution,
            pixel_errors[image_index],
            image_index
        );
#ifdef TRACK_PROGRESS
    std::cerr << "Created fit images" << std::endl;
#endif

    FitPSF::Config *fit_configuration =
        reinterpret_cast<FitPSF::Config*>(configuration);

    const PSF::Grid& grid = (
        (*fit_configuration)["psf.bicubic.grid"].as<PSF::Grid>()
    );
#ifdef TRACK_PROGRESS
    std::cerr << "Got PSF grid." << std::endl;
#endif
    PSF::PiecewiseBicubic psf(grid.x_grid.begin(),
                              grid.x_grid.end(),
                              grid.y_grid.begin(),
                              grid.y_grid.end());
#ifdef TRACK_PROGRESS
    std::cerr << "Created a PSF" << std::endl;
#endif
    std::vector<double> zeros(grid.x_grid.size() * grid.y_grid.size(),
                              0);
    psf.set_values(zeros.begin(), zeros.begin(),
                   zeros.begin(), zeros.begin());

    FitPSF::LinearSourceList fit_sources, dropped_sources;

#ifdef TRACK_PROGRESS
    std::cerr << "Created source lists" << std::endl;
#endif
    IO::H5IODataTree *real_output_data_tree =
        reinterpret_cast<IO::H5IODataTree*>(output_data_tree);

#ifdef TRACK_PROGRESS
    std::cerr << "Converted output data tree." << std::endl;
#endif
    prepare_fit_sources(
        *fit_configuration,
        fit_images,
        source_ids,
        source_coordinates,
        psf_terms,
        enabled,
        number_sources,
        number_terms,
        backgrounds,
        fit_sources,
        dropped_sources,
        subpix_map,
        psf,
        *real_output_data_tree
    );

#ifdef TRACK_PROGRESS
    std::cerr << "Got " << fit_sources.size() << " fit sources:" << std::endl;
#ifndef NDEBUG
    for(
        FitPSF::LinearSourceList::const_iterator src_i = fit_sources.begin();
        src_i != fit_sources.end();
        ++src_i
    )
        std::cerr << "x=" << (*src_i)->x() << ", y=" << (*src_i)->y()
                  << std::endl;
#endif
#endif

    Eigen::VectorXd best_fit_coef;
    FitPSF::LinearSourceList empty_source_list;

    bool ignore_dropped = (
        fit_configuration->count("psf.ignore-dropped") != 0
        &&
        (*fit_configuration)["psf.ignore-dropped"].as<bool>()
    );
#ifdef TRACK_PROGRESS
    std::cerr << "Ignore dropped: " << ignore_dropped << std::endl;
#endif
    bool converged = false;
    if(grid.x_grid.size() > 2 && grid.y_grid.size() > 2)
        converged = FitPSF::fit_piecewise_bicubic_psf(
            fit_sources,
            (ignore_dropped ? empty_source_list : dropped_sources),
            (*fit_configuration)["gain"].as<double>(),
            grid.x_grid,
            grid.y_grid,
            subpix_map,
            (*fit_configuration)[
            "psf.bicubic.max-abs-amplitude-change"
            ].as<double>(),
            (*fit_configuration)[
            "psf.bicubic.max-rel-amplitude-change"
            ].as<double>(),
            (*fit_configuration)["psf.max-chi2"].as<double>(),
            (*fit_configuration)["psf.bicubic.pixrej"].as<double>(),
            (*fit_configuration)["psf.min-convergence-rate"].as<double>(),
            (*fit_configuration)["psf.max-iterations"].as<int>(),
            (*fit_configuration)["psf.bicubic.smoothing"].as<double>(),
            best_fit_coef
        );
    else {
        converged = true;
        for(
            FitPSF::LinearSourceList::const_iterator
                src_i = fit_sources.begin();
            src_i != fit_sources.end();
            ++src_i
        ) {
            (*src_i)->flux(0).value() = 1.0;
            (*src_i)->flux(0).error() = 0.0;
            (*src_i)->flux(0).flag() = Core::GOOD;
            (*src_i)->chi2() = Core::NaN;
        }
        for(
            FitPSF::LinearSourceList::const_iterator
                src_i = dropped_sources.begin();
            src_i != dropped_sources.end();
            ++src_i
        ) {
            (*src_i)->flux(0).value() = Core::NaN;
            (*src_i)->flux(0).error() = Core::NaN;
            (*src_i)->flux(0).flag() = Core::BAD;
            (*src_i)->chi2() = Core::NaN;
        }
    }
#ifdef TRACK_PROGRESS
    std::cerr << "Converged: " << converged << std::endl;
#endif

    fit_sources.splice(fit_sources.end(), dropped_sources);
    fit_sources.sort(
        FitPSF::compare_source_assignment_ids<FitPSF::LinearSource>
    );

#ifdef TRACK_PROGRESS
    std::cerr << "Re-integrated dropped sources to source list." << std::endl;
#endif

    FitPSF::fill_output_data_tree_common(
        fit_sources,
        *real_output_data_tree,
        (*fit_configuration)["magnitude-1adu"].as<double>()
    );
#ifdef TRACK_PROGRESS
    std::cerr << "Updated common to all PSF fits entries of the output tree."
              << std::endl;
#endif
    real_output_data_tree->put("psffit.psfmap",
                               best_fit_coef,
                               IO::TranslateToAny< Eigen::VectorXd >());
#ifdef TRACK_PROGRESS
    std::cerr << "Added PSF map to output tree." << std::endl;
#endif

#ifdef VERBOSE_DEBUG
    std::cerr << "Output data tree: " << *real_output_data_tree << std::endl;
#endif

    for(
        FitPSF::LinearSourceList::iterator src_i = fit_sources.begin();
        src_i != fit_sources.end();
        ++src_i
    )
        delete *src_i;

    return converged;
}
