/**\file
 *
 * \brief Declares the base class for all sub-pixel maps.
 *
 * \ingroup Core
 */

#ifndef __SUBPIXEL_MAP_H
#define __SUBPIXEL_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Error.h"
#include <valarray>
#include <iostream>

namespace Core {

    ///\brief The base class for all sub-pixel maps.
    ///
    ///\ingroup Core
    class LIB_PUBLIC SubPixelMap {
    private:
        ///\brief The sensitivities of the sub-pixels making up the map,
        ///arranged with the x varying fast and y slow.
        std::valarray<double> __sensitivities;

        unsigned
            ///See x_resolution()
            __x_res,

            ///See x_resolution()
            __y_res;

        ///See name()
        std::string __name;
    public:
        ///\brief Create a sub-pixel map (must call set_resolution and actually
        ///set sensitivities before use).
        SubPixelMap(
            ///See name()
            const std::string &name
        ) :
            __name(name)
        {}

        ///\brief Create a sub-pixel map with the given resolution (must set
        ///sensitivities before use).
        SubPixelMap(
            ///Resolution along x
            unsigned x_res,

            ///Resolution along y
            unsigned y_res,

            ///See name()
            const std::string &name
        ) :
            __sensitivities(x_res * y_res),
            __x_res(x_res),
            __y_res(y_res),
            __name(name)
        {}

        ///Create a sub-pixel map with the given resolution and sensitivities.
        SubPixelMap(
            ///The sensitivities to set arranged with x varying faster than y.
            double *sensitivities,

            ///Resolution along x
            unsigned x_res,

            ///Resolution along y
            unsigned y_res
        ) :
            __sensitivities(sensitivities, x_res * y_res),
            __x_res(x_res),
            __y_res(y_res)
        {}

        ///Sets the resolution of the sub-pixel map, losing previous content.
        void set_resolution(unsigned x_res, unsigned y_res)
        {__x_res=x_res; __y_res=y_res; __sensitivities.resize(x_res*y_res);}

        ///Returns the resolution of the map along the x direction.
        unsigned long x_resolution() const {return __x_res;}

        ///Returns the resolution of the map along the x direction.
        unsigned long y_resolution() const {return __y_res;}

        ///\brief Returns a modifiable reference to the sensitivity at the
        ///given sub-pixel.
        double &operator()(unsigned long x, unsigned long y)
        {return __sensitivities[x+y*__x_res];}

        ///\brief Returns a copy of the sensitivity at the given sub-pixel.
        double operator()(unsigned long x, unsigned long y) const
        {return __sensitivities[x+y*__x_res];}

        ///\brief Returns the name of the given sub-pixel map (as defined at
        ///construction).
        const std::string &name() {return __name;}

        ///Copy RHS to *this.
        SubPixelMap &operator=(const SubPixelMap &RHS)
        {__sensitivities.resize(RHS.__sensitivities.size());
            __sensitivities=RHS.__sensitivities;
            __x_res=RHS.__x_res; __y_res=RHS.__y_res; return *this;}

        ///The minimum sensitivity in any part of the map.
        double min() const {return __sensitivities.min();}

        ///The maximum sensitivity in any part of the map.
        double max() const {return __sensitivities.max();}
    }; //End SubPixelMap class.

}//End Core namespace.

///Outputs the sensitivities of the subpixels as an array to the given
///stream.
std::ostream &operator<<(
    ///The stream to write to.
    std::ostream &os,

    ///The map to output.
    const Core::SubPixelMap &subpix_map
);

#endif
