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
#ifdef DEBUG
    std::cerr << "craeted config at: " << result << std::endl;
#endif
    return result;
}

void destroy_psffit_configuration(
    FittingConfiguration *configuration
)
{
    std::cerr << "destroying config at: " << configuration << std::endl;
    delete reinterpret_cast<FitPSF::Config*>(configuration);
}

void update_psffit_configuration(bool prffit,
                                 FittingConfiguration *target_configuration,
                                 ...)
{
#ifdef DEBUG
    std::cerr << "updating config at: " << target_configuration << std::endl;
#endif

    FitPSF::Config *configuration =
        reinterpret_cast<FitPSF::Config*>(target_configuration);

    char fit_mode[] = "fitpsf",
         config_file_option[] = "--config-file",
         empty_str[] = "";
    if(prffit) fit_mode[4]='r';
    char **argv_like = new char*[3];
    argv_like[0] = fit_mode;
    argv_like[1] = config_file_option;
    argv_like[2] = empty_str;

#ifdef DEBUG
    std::cerr << "Pretending executable is: " << argv_like[0] << std::endl;
#endif

    configuration->parse(3, argv_like);
    delete[] argv_like;

    opt::options_description config_file_options;

    config_file_options.add(configuration->cmdline_config_options())
                       .add(configuration->hidden_options());

    va_list arg_list;
    va_start(arg_list, target_configuration);

    for(
        char *param_name = va_arg(arg_list, char*);
        param_name[0] != '\0';
        param_name = va_arg(arg_list, char*)
    ) {
        char *param_value = va_arg(arg_list, char*);

        std::stringstream config_stream;
        config_stream << param_name << " = " << param_value << std::endl;
#ifdef DEBUG
        config_stream.seekg(0, std::ios::beg);
        std::cerr << "Setting " << config_stream.str() << std::endl;
#endif
        config_stream.seekg(0, std::ios::beg);
        opt::store(
            opt::parse_config_file(config_stream,
                                   config_file_options,
                                   true),
            *configuration
        );
    }

    va_end(arg_list);

    opt::notify(*configuration);
}

///Create a list of all fitting sources (assign pixels, filter etc.).
void prepare_fit_sources(
    ///The configuration for PSF fitting.
    const FitPSF::Config &configuration,

    ///The PSF fitting images for simultaneous fits.
    std::vector< FitPSF::Image<FitPSF::LinearSource> > fit_images,

    ///See same name argument to piecewise_bicubic_fit()
    char **column_names,

    ///See same name argument to piecewise_bicubic_fit()
    char ***source_ids,

    ///See same name argument to piecewise_bicubic_fit()
    double **column_data,

    ///See same name argument to piecewise_bicubic_fit()
    unsigned long *number_sources,

    ///See same name argument to piecewise_bicubic_fit()
    unsigned long number_columns,

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

    ///The object to add result data to (e.g. PSF map variables).
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
        image_index_stream.width(
            int(std::floor(std::log10(fit_images.size())) + 1)
        );
        image_index_stream.fill('0');
        image_index_stream << image_index;

        FitPSF::IOSources image_sources(image_index_stream.str().c_str(),
                                        source_ids[image_index],
                                        column_data[image_index],
                                        column_names,
                                        number_sources[image_index],
                                        number_columns);
#ifdef TRACK_PROGRESS
        std::cerr << "List contains "
                  << image_sources.locations().size()
                  << " sources"
                  << std::endl;
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
        FitPSF::add_expansion_terms(
            image_sources,
            configuration["psf.terms"].as<std::string>(),
            section_fit_sources,
            section_dropped_sources
        );
#ifdef TRACK_PROGRESS
        std::cerr << "Added expansion terms." << std::endl;
#endif

        if(configuration["psf.terms"].as<std::string>() != "") {
            typedef IO::IOTreeBase::path_type path;
            output_data_tree.put(
                path(
                    "psffit|variables|" + image_index_stream.str(),
                    '|'
                ),
                image_sources.columns(),
                IO::TranslateToAny<PSF::MapVarListType>()
            );
        }
#ifdef TRACK_PROGRESS
        std::cerr << "Added PSF fit variables to result tree." << std::endl;
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

bool piecewise_bicubic_fit(double **pixel_values,
                           double **pixel_errors,
                           char **pixel_masks,
                           unsigned long number_images,
                           unsigned long image_x_resolution,
                           unsigned long image_y_resolution,
                           char **column_names,
                           char ***source_ids,
                           double **column_data,
                           unsigned long *number_sources,
                           unsigned long number_columns,
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
            pixel_errors[image_index]
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
        column_names,
        source_ids,
        column_data,
        number_sources,
        number_columns,
        backgrounds,
        fit_sources,
        dropped_sources,
        subpix_map,
        psf,
        *real_output_data_tree
    );

#ifdef TRACK_PROGRESS
    std::cerr << "Got " << fit_sources.size() << " fit sources:" << std::endl;
    for(
        FitPSF::LinearSourceList::const_iterator src_i = fit_sources.begin();
        src_i != fit_sources.end();
        ++src_i
    )
        std::cerr << "x=" << (*src_i)->x() << ", y=" << (*src_i)->y()
                  << std::endl;
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
    bool converged = FitPSF::fit_piecewise_bicubic_psf(
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

    return converged;
}
