///\ingroup SubPixPhot
///\ingroup FitSubpix

/**\file
 *
 * \brief Declares the base class for background extract(ion/ed)
 * from an astronomical image.
 *
 * \ingroup Background
 */

#ifndef __SOURCE_BACKGROUND_H
#define __SOURCE_BACKGROUND_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Image.h"
#include "../Core/Point.h"
#include "../Core/NaN.h"
#include <list>
#include <limits>

namespace Background {

    ///The base class for the backgnound under a source.
    class LIB_LOCAL Source {
    private:
        double __value, ///< An estimate of the background for the source
               __error; ///< An estimate of the uncertainty in the value.
        unsigned __pixels; ///< How many pixels were incuded in the determination
    public:
        ///Create a background object.
        Source(
            ///Value of the background
            double value = Core::NaN,

            ///Uncertainty in the background value
            double error = Core::NaN,

            ///Number of pixels used in the measurement
            unsigned pixels=0
        ) :
            __value(value), __error(error), __pixels(pixels) {}


        ///@{
        ///\brief The best estimate of the background level under the source.
        double value() const {return __value;}
        double &value() {return __value;}
        ///@}


        ///@{
        ///The best estimate of the uncertainty in the backgronud determination.
        double error() const {return __error;}
        double &error() {return __error;}
        ///@}

        ///@{
        ///The number of pixels that contributed to the determination of the
        ///value and error.
        unsigned pixels() const {return __pixels;}
        unsigned &pixels() {return __pixels;}
        ///@}
    }; //End Source class.

    ///Output Human readable representation of the source background to stream.
    std::ostream &operator<<(
        ///The stream to write to.
        std::ostream &os,

        ///The background under a source to output.
        const Source& source
    );

}//End Core namespace.

#endif

