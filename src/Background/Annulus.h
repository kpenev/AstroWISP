/**\file
 *
 * \brief Declare a class defining an annulus around each source to measure
 * the background in.
 *
 * \ingroup Background
 */

#ifndef __BACKGROUND_ANNULUS
#define __BACKGROUND_ANNULUS

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/NaN.h"
#include <string>
#include <sstream>

namespace Background {

    ///\brief The annulus used to measure the background under a source.
    class LIB_PUBLIC Annulus {
    private:
        double __inner_radius, ///< The radius of the inside boundary
               __width; ///< Outer radius=width+inner radius
    public:
        ///Make an annulus with the given inner radius and width.
        Annulus(
            double inner_radius = Core::NaN,
            double width = Core::NaN
        ) :
            __inner_radius(inner_radius),
            __width(width)
        {}

        ///The radius of the inside boundary
        double inner_radius() const {return __inner_radius;}

        ///The radius of the inside boundary
        double &inner_radius() {return __inner_radius;}
        
        ///The radius of the outside boundary
        double outer_radius() const {return __inner_radius+__width;}

        ///The difference between the radi of the inside and outside boundaries
        double width() const {return __width;}

        ///The difference between the radi of the inside and outside boundaries
        double &width() {return __width;}

        ///Represent the annulus as a string: '\<inner\>:\<width\>'.
        operator std::string() const;
    }; //End Annulus class.

} //End Background namespace.

#endif
