/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of
 * SubPixPhot.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SUBPIX_PHOT_C_INTERFACE_H
#define __SUBPIX_PHOT_C_INTERFACE_H

#include "../IO/CInterface.h"
#include "../Core/CInterface.h"

extern "C" {
    ///Opaque struct to cast to/from SubPixPhot::Config.
    struct LIB_PUBLIC SubPixPhotConfiguration;

    ///\brief Create an object for holding the configuration on how to perform
    ///aperture photometry.
    //
    ///Must be initialized by updated_subpixphot_configuration() before use.
    ///
    ///Must be destroyed by destroy_subpixphot_configuration when no longer needed
    ///to avoid a memory leak.
    LIB_PUBLIC SubPixPhotConfiguration *create_subpixphot_configuration();

    ///\brief Destroy a configuration previously created by
    ///create_subpixphot_configuration()
    LIB_PUBLIC void destroy_subpixphot_configuration(
        ///The configuration to destroy.
        SubPixPhotConfiguration *configuration
    );

    ///Update the configuration for aperture photometry.
    LIB_PUBLIC void update_subpixphot_configuration(
        ///The configuration to update.
        SubPixPhotConfiguration *target_configuration,

        ///Alternating <parameter name>, <parameter value> pairs, with both
        ///etries being of type char* type.
        ...
    );

    ///\brief Measure the flux of sources in an image using sub-pixel aware
    ///aperture photometry.
    LIB_PUBLIC void subpixphot(
        ///The image to perform aperture photometry on.
        const CoreImage *image,

        ///The sub-pixel map to assume.
        const CoreSubPixelMap *subpixmap,

        ///The configuration for how to performe the aperture photometry.
        SubPixPhotConfiguration *configuration,

        ///The tree containing the resuts from PSF fitting. On exit, updated
        ///with the aperture photometry results.
        H5IODataTree *io_data_tree,

        ///The image index within the output tree image corresponds to.
        unsigned image_index
    );

} //End Extern "C"

#endif
