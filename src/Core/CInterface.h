/**\file
 *
 * \brief Declare some opaque structures to be used when creating C-interfaces
 * to libraries.
 *
 * \ingroup Core
 */

extern "C" {

    ///Opaque struct to cast to/from Core::Image.
    struct CoreImage;

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
};//End extern "C".
