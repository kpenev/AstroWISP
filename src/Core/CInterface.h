/**\file
 *
 * \brief Declare some opaque structures to be used when creating C-interfaces
 * to libraries.
 *
 * \ingroup Core
 */

#ifndef __CORE_C_INTERFACE_H
#define __CORE_C_INTERFACE_H

#include "SharedLibraryExportMacros.h"

extern "C" {
    ///Opaque struct to cast to/from Core::Image.
    struct LIB_PUBLIC CoreImage;

    ///Opaque struct to cast to/from Core::Flux.
    struct CoreFlux;

    ///Opaque struct to cast to/from Core::FluxPair.
    struct CoreFluxPair;

    ///Opaque struct to cast to/from Core::Point.
    struct CorePoint;

    ///Opaque struct to cast to/from Core::SourceID.
    struct CoreSourceID;

    ///Opaque struct to cast to/from Core::SourceLocation.
    struct CoreSourceLocation;

    ///Opaque struct to cast to/from Core::SubPixelCorrectedFlux.
    struct CoreSubPixelCorrectedFlux;

    ///Opaque struct to cast to/from Core::SubPixelMap.
    struct CoreSubPixelMap;

    ///Create and fill image for processing by AstroWISP's tools.
    LIB_PUBLIC CoreImage *create_core_image(
        ///See same name argument to Core::Image::Image().
        unsigned long x_resolution,

        ///See same name argument to Core::Image::Image().
        unsigned long y_resolution,

        ///See same name argument to Core::Image::Image().
        double *values,

        ///See same name argument to Core::Image::Image().
        double *errors,

        ///See same name argument to Core::Image::Image(). May be NULL to
        ///disable pixel masking.
        char *mask,

        ///Should the image simply wrap the given data (makes a copy if false).
        ///If true, the input data should not be released bofore this object is
        ///destroyed.
        bool wrap
    );

    ///Relese memory for an image created by create_core_image().
    LIB_PUBLIC void destroy_core_image(
        ///The image to destroy. Must have been created by create_core_image().
        CoreImage *image
    );

    ///Create and fill a sub-pixel sensitivity map.
    LIB_PUBLIC CoreSubPixelMap *create_core_subpixel_map(
        ///See x_res argument to Core::SubPixelMap::SubPixelMap().
        unsigned long x_resolution,

        ///See y_res argument to Core::SubPixelMap::SubPixelMap().
        unsigned long y_resolution,

        ///The sensitivities for each piece of the map. The first x_resolution
        ///entries give the sensitivity of pixels at y=0 for x going from 0 to
        ///x_resolution. Followed by the sensitivites at y=1 etc.
        double *sensitivities
    );

    ///Release memory for a sub-pixel map created by create_core_subpixel_map().
    LIB_PUBLIC void destroy_core_subpixel_map(
        ///The sub-pixel map to destroy. Must have been created by
        ///create_core_subpixel_map()
        CoreSubPixelMap *map
    );
};//End extern "C".

#endif
