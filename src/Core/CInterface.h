/**\file
 *
 * \brief Declare some opaque structures to be used when creating C-interfaces
 * to libraries.
 *
 * \ingroup Core
 */

#include "SharedLibraryExportMacros.h"
#include "Image.h"

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
    LIB_PUBLIC void destroy_core_image(CoreImage *image);
};//End extern "C".
