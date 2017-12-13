/**\file
 *
 * \brief Define the backrgound excess functions and some methods of Pixel.h
 *
 * \ingroup FitPSF
 */

#include "Pixel.h"
#include "Source.h"

namespace FitPSF {

    double background_excess(double value,
                             double variance,
                             double background_value,
                             double background_variance)
    {
        return (std::max(value - background_value, 0.0)
                /
                std::sqrt(variance + background_variance));
    }

    double background_excess(double                     value,
                             double                     variance,
                             const Background::Source   &background,
                             double                     gain)
    {
        return background_excess(value,
                                 variance,
                                 background.value() * gain,
                                 std::pow(background.error() * gain, 2));
    }

} //End FitPSF namespace.
