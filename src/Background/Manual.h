/**\file
 *
 * \brief Declare a background class for sources with known backgrounds.
 *
 * \ingroup Background
 */

#ifndef __MANUAL_BACKGROUND_H
#define __MANUAL_BACKGROUND_H

#include "Measure.h"

namespace Background {

    ///\brief The background is manually specified for a list of sources.
    class LIB_LOCAL Manual : public Measure {
    private:
        ///The sources for which background is known.
        std::vector<Source> __sources;

        ///The current source.
        std::vector<Source>::const_iterator __current_source;
    public:
        ///Construct from c-style arrays.
        Manual(
            ///The background values for all sources.
            double *value,

            ///Estimate of the background error for all sources.
            double *error,

            ///How many pixels is each of the value and error entries based on.
            unsigned *num_pixels,

            ///How many sources are known (the length of all of the above
            ///arrays).
            unsigned num_sources
        );


        ///Adding sources is not supported.
        virtual void add_source(double, double)
        {
            throw Error::NotImplemented(
                "Adding sources to Backgronud::Manual is not supported!"
            );
        }

        ///Adding sources is not supported.
        virtual void add_source(const Core::Point<double> &)
        {
            throw Error::NotImplemented(
                "Adding sources to Backgronud::Manual is not supported!"
            );
        }

        ///Estimate the background around the current source.
        virtual Source operator()() const
        {return *__current_source;}

        ///\brief Jump to the next source.
        ///
        ///The next evaluation (operator()) will return the background 
        ///estimate for the next source. Retruns true iff the new source
        ///exists.
        virtual bool next_source()
        {return ++__current_source == __sources.end();}

        ///Move one source back. Returns true iff the new source exists.
        virtual bool previous_source()
        {
            if(__current_source == __sources.begin()) return false;
            else {
                --__current_source;
                return true;
            }
        }

        ///Restart iteration over sources from the beginning.
        virtual void jump_to_first_source()
        {__current_source = __sources.begin();}

    };//End Manual class.
} //End Background namespace

#endif
