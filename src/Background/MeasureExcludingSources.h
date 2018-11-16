/**\file
 *
 * \brief Declares a class for extracting the background around sources
 * attempting to exclude nearby sources.
 *
 * \ingroup Background
 */

#ifndef __MEASURE_EXCLUDING_SOURCES_H
#define __MEASURE_EXCLUDING_SOURCES_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Measure.h"
#include "../Core/Image.h"
#include "../Core/Point.h"
#include "../Core/NaN.h"
#include <list>
#include <limits>
#include <cmath>

namespace Background {

    ///\brief Avoid pixels close to sources when extracting a background.
    ///
    ///A base class for background extractors which exclude some circular
    ///aperture around every source when determining the backgrounds of other
    ///sources.
    ///
    ///Masks (assigns NaN to) the pixels in the given image whose centers lie
    ///within the specified aperture of at least one source. This masked array
    ///can then be used by a non-abstract child class to derive a background
    ///value and uncertainty.
    class LIB_PUBLIC MeasureExcludingSources : public Measure {
    private:
        ///The radius of the exclusion area around each source.
        double __aperture;

        ///An iterator over the sources in the image.
        std::list< Core::Point<double> >::const_iterator __source_iter;

        ///Excludes all pixels which have a mask other than GOOD.
        void exclude_mask();

    protected:
        ///The sources in the image
        std::list< Core::Point<double> > _sources;

        ///The original image with the areas near sources masked.
        Core::Image<double> _bg_values;
    public:
        ///\brief Construct a background extractor for the given image excluding
        ///pixels with centers inside the given aperture around all sources.
        MeasureExcludingSources(
            ///The image for which background extraction is to be done.
            const Core::Image<double> &image,

            ///The size of the aperture to exclude around each source.
            double exclude_aperture
        ) :
            __aperture(exclude_aperture),
            _bg_values(image)
        {
            exclude_mask();
        }

        ///\brief Construct a background extractor for the given image excluding
        ///pixels with centers inside the given aperture around all listed
        ///sources.
        template<class POINT_TYPE>
        MeasureExcludingSources(
            ///The image for which background extraction is to be done.
            const Core::Image<double> &image,

            ///The size of the aperture to exclude around each source.
            double exclude_aperture,

            ///The sources in the image.
            const std::list< POINT_TYPE > &sources
        );

        ///\brief Add a source.
        ///
        ///Register another source and remove the values in its aperture from
        ///_bg_values.
        virtual void add_source(double x, double y);

        ///See add_source(double, double)
        virtual void add_source(const Core::Point<double> &location)
        {add_source(location.x(), location.y());}

        ///See add_source(double, double)
        virtual void add_source(const Core::Point<double> *location)
        {add_source(location->x(), location->y());}

        ///Estimate the background around the current source.
        inline virtual Source operator()() const
        {return operator()(*__source_iter);}

        ///Estimate the background around the given position.
        virtual Source operator()(const Core::Point<double> &location) const
        {return operator()(location.x(), location.y());}

        ///Estimate the background around the given position.
        virtual Source operator()(double x, double y) const =0;

        ///\brief Jump to the next source and return true iff it exists.
        ///
        ///The next evaluation (operator()) will return the background
        ///estimate for the next source. Retruns true iff the new source exists.
        inline virtual bool next_source()
        {return ++__source_iter!=_sources.end();}

        ///Move one source back. Returns true iff the new source exists.
        inline virtual bool previous_source()
        {
            if(__source_iter==_sources.begin()) return false;
            else {--__source_iter; return true;}
        }

        ///Restart iteration over sources from the beginning.
        virtual void jump_to_first_source()
        {__source_iter=_sources.begin();}
    }; //End MeasureExcludingSources class.

    template<class POINT_TYPE>
    MeasureExcludingSources::MeasureExcludingSources(
        const Core::Image<double> &image,
        double exclude_aperture,
        const std::list< POINT_TYPE > &sources
    ) :
        __aperture(exclude_aperture),
        _bg_values(image)
    {
        exclude_mask();
        for(
            typename std::list< POINT_TYPE >::const_iterator
            src_iter=sources.begin(); src_iter!=sources.end();
            src_iter++
        )
            add_source(*src_iter);
    }

}//End Background namespace.

#endif
