/**\file
 *
 * \brief Declare a class for flux measurements.
 *
 * \ingroup Core
 */

#ifndef __FLUX_H
#define __FLUX_H

#include "Typedefs.h"

namespace Core {

    ///\brief A class representing the flux measurement for a source, 
    ///including an error estimate and a flag.
    class Flux {
    private:
        double __value, ///< The estimate of the flux.
               __error; ///< An estimate of the error.
        Core::PhotometryFlag __flag; ///< The quality flag.
    public:
        Flux(
            double value = NaN,
            double error = NaN,
            Core::PhotometryFlag flag = UNDEFINED
        ) :
            __value(value),
            __error(error),
            __flag(flag)
        {}
        double value() const {return __value;}
        double &value() {return __value;}
        double error() const {return __error;}
        double &error() {return __error;}
        Core::PhotometryFlag flag() const {return __flag;}
        Core::PhotometryFlag &flag() {return __flag;}
    };

} //End Core namespace.
#endif
