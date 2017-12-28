/**\file
 *
 * \brief Defines a point class in two dimensions.
 *
 * \ingroup Core
 */

#ifndef __POINT_H
#define __POINT_H

#include "../Core/SharedLibraryExportMacros.h"

namespace Core {

    ///\brief Representing a 2D point with real valued cartesian coordinates.
    template <class COORD_TYPE>
        class LIB_PUBLIC Point {
        private:
            COORD_TYPE
                ///The abscissa coordinate of the point.
                __x,

                ///The oordinate coordinate of the point.
                __y;
        public:
            ///Create a point at the given location.
            Point(COORD_TYPE x=0, COORD_TYPE y=0) : __x(x), __y(y) {}

            ///The abscissa coordinate of the point.
            COORD_TYPE x() const {return __x;}

            ///The abscissa coordinate of the point.
            COORD_TYPE &x() {return __x;}

            ///The oordinate coordinate of the point.
            COORD_TYPE y() const {return __y;}

            ///The oordinate coordinate of the point.
            COORD_TYPE &y() {return __y;}
        }; //End Point class.

} //End Core namespace.

#endif
