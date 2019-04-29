#include "CInterface.h"
#include "Config.h"
#include "Common.h"
#include "../PSF/PiecewiseBicubicMap.h"
#include "../Core/SubPixelMap.h"
#include "../Core/Image.h"
#include "../Core/SubPixelCorrectedFlux.h"
#include "../PSF/DataTreeCalculations.h"

#include <cstdarg>
#include <iostream>

SubPixPhotConfiguration *create_subpixphot_configuration()
{
    SubPixPhotConfiguration* result =
        reinterpret_cast<SubPixPhotConfiguration*>(
            new SubPixPhot::Config(0, NULL)
        );
#ifdef VERBOSE_DEBUG
    std::cerr << "craeted subpixphot config at: " << result << std::endl;
#endif
    return result;
}

void destroy_subpixphot_configuration(SubPixPhotConfiguration *configuration)
{
#ifdef VERBOSE_DEBUG
    std::cerr << "destroying subpixphot config at: "
              << configuration
              << std::endl;
#endif
    delete reinterpret_cast<SubPixPhot::Config*>(configuration);
}

void update_subpixphot_configuration(
    SubPixPhotConfiguration *target_configuration,
    ...
)
{
    SubPixPhot::Config *configuration =
        reinterpret_cast<SubPixPhot::Config*>(target_configuration);

    va_list options;
    va_start(options, target_configuration);

    IO::update_configuration(*configuration, "subpixphot", options);

    va_end(options);
}

LIB_PUBLIC void subpixphot(const CoreImage *image,
                           const CoreSubPixelMap *subpixmap,
                           SubPixPhotConfiguration *configuration,
                           H5IODataTree *io_data_tree,
                           unsigned image_index)
{
#ifdef TRACK_PROGRESS
    std::cerr << "Starting aperture photometry." << std::endl;
#endif
    const Core::SubPixelMap *real_subpixmap =
        reinterpret_cast<const Core::SubPixelMap *>(subpixmap);

    const Core::Image<double> *real_image =
        reinterpret_cast<const Core::Image<double> *>(image);

    const SubPixPhot::Config *real_configuration =
        reinterpret_cast<const SubPixPhot::Config *>(configuration);

    std::ostringstream image_index_str;
    image_index_str << image_index;

    IO::H5IODataTree *real_io_data_tree =
        reinterpret_cast<IO::H5IODataTree*>(io_data_tree);

    PSF::fill_psf_amplitudes(*real_io_data_tree, image_index_str.str());

    std::cerr << "Setting apertures" << std::endl;
    Core::RealList apertures =
        (*real_configuration)["ap.aperture"].as<Core::RealList>();
    apertures.sort();

    std::cerr << "Reading PSF map" << std::endl;
    PSF::PiecewiseBicubicMap psf_map(*real_io_data_tree,
                                     apertures.back() + 1.0);

    std::cerr << "Creating flux measuring object" << std::endl;
    Core::SubPixelCorrectedFlux<Core::SubPixelMap> measure_flux(
        *real_image,
        *real_subpixmap,
        (*real_configuration)["ap.const-error"].as<double>(),
        apertures,
        (*real_configuration)["gain"].as<double>()
    );

    std::cerr << "Measuring flux for image "
              << image_index_str.str()
              << std::endl;
    SubPixPhot::add_flux_measurements(
        psf_map,
        measure_flux,
        (*real_configuration)["magnitude-1adu"].as<double>(),
        *real_io_data_tree,
        image_index_str.str()
    );
}
