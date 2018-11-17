/**\file
 *
 * \brief Declare a class for flux measurements.
 *
 * \ingroup Core
 */

#ifndef __FLUX_H
#define __FLUX_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Typedefs.h"

namespace Core {

    ///\brief A class representing the flux measurement for a source, 
    ///including an error estimate and a flag.
    class LIB_PUBLIC Flux {
    private:
        double __value, ///< See value()
               __error; ///< See error()
        Core::PhotometryFlag __flag; ///< See flag()
    public:
        ///Create a flux with the given properties.
        Flux(
            ///The flux value
            double value = NaN,

            ///Error estimate for the flux value.
            double error = NaN,

            ///Quality flag indicating the relaibility of the flux measurement.
            Core::PhotometryFlag flag = UNDEFINED
        ) :
            __value(value),
            __error(error),
            __flag(flag)
        {}

        ///The estimate of the flux (immutable).
        double value() const {return __value;}

        ///The estimate of the flux (mutable).
        double &value() {return __value;}

        ///An estimate of the error (immutable).
        double error() const {return __error;}

        ///An estimate of the error (mutable).
        double &error() {return __error;}

        ///The quality flag (immutable)
        Core::PhotometryFlag flag() const {return __flag;}

        ///The quality flag (mutable)
        Core::PhotometryFlag &flag() {return __flag;}
    };

} //End Core namespace.
#endif
