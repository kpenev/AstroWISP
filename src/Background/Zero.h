/**\file
 * 
 * \brief Declare a background class for images from which the background has
 * already been subtracted.
 *
 * \ingroup Background
 */

#ifndef __ZERO_BACKGROUND_H
#define __ZERO_BACKGROUND_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Point.h"
#include "Measure.h"

namespace Background {

    ///\brief The background is assumed already subtracted from the input 
    ///image.
    class LIB_LOCAL Zero : public Measure {
    private:
        ///How many sources have been added so far.
        unsigned __number_sources;

        ///The index of the currently selected source.
        unsigned __current_index;
    public:
        ///Default constructor.
        Zero() : __number_sources(0), __current_index(0) {}

        ///Zero background for the given number of sources.
        Zero(unsigned number_sources) :
            __number_sources(number_sources),
            __current_index(0)
        {}

        ///See Measure::add_source(double, double)
        void add_source(double, double) {++__number_sources;}

        ///See Measure::add_source(const Core::Point<double> &)
        void add_source(const Core::Point<double> &) {++__number_sources;}

        ///Estimate the background around the current source.
        Source operator()() const 
        {
            return Source(0.0, 0.0, std::numeric_limits<unsigned>::max());
        }

        ///\brief Jump to the next source.
        ///
        ///The next evaluation (operator()) will return the background 
        ///estimate for the next source. Retruns true iff the new source 
        ///exists.
        bool next_source() 
        {
            assert(__current_index < __number_sources);
            ++__current_index;
            return __current_index < __number_sources;
        }

        ///Move one source back. Returns true iff the new source exists.
        bool previous_source()
        {
            if(__current_index > 0) {--__current_index; return true;}
            else return false;
        }

        ///Restart iteration over sources from the beginning.
        virtual void jump_to_first_source() {__current_index = 0;}
    }; //End Zero class.

} //End Background namespace.

#endif
