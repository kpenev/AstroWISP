/**\file
 *
 * \brief Declares the base class for all sub-pixel maps.
 *
 * \ingroup Core
 */

#ifndef __SUBPIXEL_MAP_H
#define __SUBPIXEL_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "fitsio.h"
#include "Error.h"
#include <valarray>
#include <iostream>

namespace Core {

    ///\brief The base class for all sub-pixel maps.
    ///
    ///\ingroup Core
    class LIB_PUBLIC SubPixelMap {
    private:
        std::valarray<double> __sensitivities;
        unsigned __x_res, __y_res;
        std::string __name;
    public:
        ///\brief Create a sub-pixel map (must call set_resolution and actually
        ///set sensitivities before use).
        SubPixelMap(const std::string &name) : __name(name) {}

        ///\brief Create a sub-pixel map with the given resolution (must set
        ///sensitivities before use).
        SubPixelMap(unsigned x_res, unsigned y_res, const std::string &name) :
            __sensitivities(x_res*y_res), __x_res(x_res), __y_res(y_res),
            __name(name) {}

        ///Create a sub-pixel map with the given resolution and sensitivities.
        SubPixelMap(double *sensitivities, unsigned x_res, unsigned y_res) :
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

        ///Saves the sub-pixel map to a FITS image with the given name.
        ///
        ///If a file with the given name already exists and the filename does 
        ///not start with '!', an exception is raised.
        void save_to_fits(const std::string &filename);

        ///The minimum sensitivity in any part of the map.
        double min() const {return __sensitivities.min();}

        ///The maximum sensitivity in any part of the map.
        double max() const {return __sensitivities.max();}
    }; //End SubPixelMap class.

}//End Core namespace.

///Outputs the sensitivities of the subpixels as an array to the given 
///stream.
std::ostream &operator<<(std::ostream &os,
                         const Core::SubPixelMap &subpix_map);

#endif
