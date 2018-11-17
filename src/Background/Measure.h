/**\file
 *
 * \brief Declares the minimum requirements for a background extractor class.
 *
 * \ingroup Background
 */

#ifndef __MEASURE_BACKGROUND_H
#define __MEASURE_BACKGROUND_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Source.h"
#include "../Core/Point.h"

namespace Background {

    ///The minimum requirements for a Background extractor class.
    class LIB_LOCAL Measure {
    public:
        ///\brief Notify the background extractor of (another) source at the
        ///given location
        virtual void add_source(double x, double y) = 0;

        ///See add_source(double, double).
        virtual void add_source(const Core::Point<double> &location) = 0;

        ///Estimate the background around the current source.
        virtual Source operator()() const =0;

        ///\brief Jump to the next source.
        ///
        ///The next evaluation (operator()) will return the background 
        ///estimate for the next source. Retruns true iff the new source
        ///exists.
        virtual bool next_source()=0;

        ///Move one source back. Returns true iff the new source exists.
        virtual bool previous_source()=0;

        ///Restart iteration over sources from the beginning.
        virtual void jump_to_first_source()=0;

        ///Do nothing by default.
        virtual ~Measure() {}
    };

} //End Background namespace.

#endif
