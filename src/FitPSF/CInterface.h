/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of
 * FitPSF.
 *
 * \ingroup FitPSF
 */

#ifndef __FITPSF_C_INTERFACE_H
#define __FITPSF_C_INTERFACE_H

#include "../Core/CInterface.h"
#include "../IO/CInterface.h"
#include "../Background/CInterface.h"

extern "C" {
    ///Opaque struct to cast to/from FitPSF::Config.
    struct LIB_PUBLIC FittingConfiguration;

    ///Opaque struct to cast to/from FitPSF::Image.
    struct LIB_PUBLIC FitPSFImage;

    ///\brief Create an object for holding the configuration for PSF fitting.
    ///
    ///Must be initialized by updated_psffit_configuration() before use.
    ///
    ///Must be destroyed by destroy_psffit_configuration when no longer needed
    ///to avoid a memory leak.
    LIB_PUBLIC FittingConfiguration *create_psffit_configuration();

    ///\brief Destroy a configuration previously created by
    ///create_psffit_configuration()
    LIB_PUBLIC void destroy_psffit_configuration(
        ///The configuration to destroy.
        FittingConfiguration *configuration
    );

    ///Update the configuration for PSF fitting.
    LIB_PUBLIC void update_psffit_configuration(
        ///Are we deing PRF fittindg (in liu of PSF fitting)?
        bool prffit,

        ///The configuration to update.
        FittingConfiguration *target_configuration,

        ///Alternating <parameter name>, <parameter value> pairs, with both
        ///etries being of type char* type.
        ...
    );

    ///Fit a smoothly varying piecewise bicubic PSF model.
    LIB_PUBLIC bool piecewise_bicubic_fit(
        ///The pixel values of the calibrated images to perform simultaneous
        ///fits of.
        double **pixel_values,

        ///Estimated pixel errors of the calibrated images to perform
        ///simultaneous fits of.
        double **pixel_errors,

        ///Mask flags of the calibrated images to perform simultaneous fits of.
        char **pixel_masks,

        ///How many images are being simultaneously fit.
        unsigned long number_images,

        ///The common x resolution of the images being processed.
        unsigned long image_x_resolution,

        ///The common y resolution of the images being processed.
        unsigned long image_y_resolution,

        ///The IDs to assign to the sources in each image. The first index is
        ///the image index the second index is the source index within each
        ///image.
        char ***source_ids,

        ///The cordinates of the source centroids. The first index is the image
        ///index, the remaining array is organized as described by the
        ///source_coordinates argument of FitPSF::IOSources::IOSources().
        double **source_coordinates,

        ///The values of each term PSF parameters depend on for each source. The
        ///first index is the image index, the remaining array is organized as
        ///described by the psf_terms argument of
        ///FitPSF::IOSources::IOSources().
        double **psf_terms,

        ///Flag indicating whether each source is allowed to participate in
        ///shape fitting or not. The first index is the image index, the
        ///remaining array is organized as described by the
        ///enabled argument of FitPSF::IOSources::IOSources().
        bool **enabled,

        ///How many sources are in column_data for each image.
        unsigned long *number_sources,

        ///How many terms are used in PSF fitting.
        unsigned long number_terms,

        ///Pointers to the measured background for the input sources indexed by
        ///the image index.
        BackgroundMeasureAnnulus **backgrounds,

        ///The configuration for how to do the fitting.
        FittingConfiguration *configuration,

        ///The sensitivities of the pixel at each sub-pixel part.
        double *subpix_sensitivities,

        ///Subpixmap x resolution.
        unsigned long subpix_x_resolution,

        ///Subpixmap y resolution.
        unsigned long subpix_y_resolution,

        ///A data tree to fill with the results.
        H5IODataTree *output_data_tree
    );

};//End extern "C"

#endif
