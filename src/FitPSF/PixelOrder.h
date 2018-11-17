/**\file
 *
 * \brief Define a class that orders pixels inside a source.
 *
 * \ingroup FitPSF
 */

#ifndef __PIXEL_ORDER_H
#define __PIXEL_ORDER_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Pixel.h"

namespace FitPSF {

    ///\brief Pixel comparison predicate.
    class LIB_LOCAL PixelOrder
    {
    private:
        double
            ///See Source::__background_electrons.
            __background_electrons,

            ///See Source::__background_electrons_variance.
            __background_electrons_variance;
    public:
        ///Construct the predicate based on a given background.
        PixelOrder(
            ///The number of electrons per pixel in the background
            double background_electrons,

            ///The estiamted variance in background_electrons.
            double background_electrons_variance
        ) :
            __background_electrons(background_electrons),
            __background_electrons_variance(background_electrons_variance)
        {}

        ///\brief Return True iff the first pixel belongs bofore the second
        ///in the ordered list of pixels.
        template<class FIT_SOURCE_TYPE>
        bool operator()(
            const Pixel<FIT_SOURCE_TYPE> *first,
            const Pixel<FIT_SOURCE_TYPE> *second
        ) const;
    };//End PixelOrder class definition.

    template<class FIT_SOURCE_TYPE>
        bool PixelOrder::operator()(
            const Pixel<FIT_SOURCE_TYPE> *first,
            const Pixel<FIT_SOURCE_TYPE> *second
        ) const
        {
            if(first->shape_fit()) {
                if(!second->shape_fit()) return true;
            } else if(first->flux_fit()) {
                if(second->shape_fit()) return false;
                else if(!second->flux_fit()) return true;
            } else if(second->shape_fit() || second->flux_fit()) return false;

            if(first->flag() < second->flag())
                return true;
            else
                return (
                    background_excess(*first,
                                      __background_electrons,
                                      __background_electrons_variance)
                    >
                    background_excess(*second,
                                      __background_electrons,
                                      __background_electrons_variance)
                );
        }

} //End FitPSF namespace.

#endif
