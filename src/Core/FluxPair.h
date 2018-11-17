/**\file
 *
 * \brief Define a class for holding subpixel sensitivity weighted and
 * non-weighted fluxes.
 *
 * \ingroup SubPixPhot
 */

#ifndef __FLUX_PAIR_H
#define __FLUX_PAIR_H

#include "../Core/SharedLibraryExportMacros.h"

namespace Core {

    ///\brief A class with two flux components corresponding to the flux
    ///corrected for subpixel sensitivity variations and partially being
    ///within the aperture and the full flux over a pixel from the assumed
    ///PSF.
    class LIB_LOCAL FluxPair {
    private:
        double
            ///See raw_flux()
            __raw_flux,

            ///See weighted_flux()
            __weighted_flux;

    public:
        ///Create a pair with the given weigthed and unweigthed fluxes.
        FluxPair(
            ///See raw_flux()
            double raw = 0,

            ///See weighted_flux()
            double weighted = 0
        ) :
            __raw_flux(raw),
            __weighted_flux(weighted)
        {}

        ///A simple integral of the PSF over the area of a 1x1 pixel.
        double raw_flux() const {return __raw_flux;}

        ///A simple integral of the PSF over the area of a 1x1 pixel.
        double &raw_flux() {return __raw_flux;}

        ///\brief An of the subpixel sensitivity times the PSF over the part
        ///of the pixel that fits within the aperture.
        double weighted_flux() const {return __weighted_flux;}

        ///\brief An of the subpixel sensitivity times the PSF over the part
        ///of the pixel that fits within the aperture.
        double &weighted_flux() {return __weighted_flux;}
    }; //End FluxPair class.

} //End SubPixPhot namespace.

#endif
