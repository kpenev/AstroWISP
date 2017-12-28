/**\file
 *
 * \brief Declares a class defining the minimum requirements for specifying
 * the locations of sources.
 *
 * \ingroup Core
 */

#ifndef __SOURCE_LOCATION_H
#define __SOURCE_LOCATION_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SourceID.h"
#include "Point.h"
#include "NaN.h"

namespace Core {

    ///The minimum requirements for objects specifying where a source is.
    class LIB_LOCAL SourceLocation : public Point<double> {
    private:
        ///An identifying string for the source.
        SourceID __id;
    public:

        ///Create a source at the given location (x,y).
        SourceLocation(double x = NaN, double y = NaN)
        : Point<double>(x, y)
        {}

        ///Create a source with the given id at the given location (x,y).
        SourceLocation(
            const SourceID& id,
            double          x = NaN,
            double          y = NaN
        )
        :
            Point<double>(x, y),
            __id( id )
        {}

        ///Identifying string for this source.
        const SourceID &id() const {return __id;}

        ///Identifying string for this source.
        SourceID &id() {return __id;}

        ///Check if the given source has the same ID as this.
        bool operator==(const SourceLocation &rhs) const
        {return __id == rhs.__id;}
    }; //End SourceLocation class.

} //End Core namespace.

#endif
