/**\file
 *
 * \brief The implementation of the methods of the MeasureExcludingSources
 * class.
 *
 * \ingroup Background
 */

#include "MeasureExcludingSources.h"
#include <iostream>

namespace Background {

    void MeasureExcludingSources::exclude_mask()
    {
        for (unsigned long yi = 0; yi < _bg_values.y_resolution(); ++yi)
            for (
                unsigned long xi = 0;
                xi < _bg_values.x_resolution();
                ++xi
            ) {
                if(!_bg_values.good(xi, yi))
                    _bg_values(xi, yi) = Core::NaN;
            }
    }

    void MeasureExcludingSources::add_source(double x, double y)
    {
        unsigned long
            xi_min = std::max(0.0, std::ceil(x - __aperture - 0.5)),
            xi_max=std::min(
                static_cast<double>(_bg_values.x_resolution()), 
                std::ceil(x+__aperture-0.5)
            ),
            yi_min=std::max(0.0, std::ceil(y-__aperture-0.5)),
            yi_max=std::min(
                static_cast<double>(_bg_values.y_resolution()), 
                std::ceil(y+__aperture-0.5)
            );
        double ap2 = __aperture * __aperture;
        _sources.push_back(Core::Point<double>(x, y));

        double xcent = x - 0.5;
        double ycent = y - 0.5;
#ifdef VERBOSE_DEBUG
        std::cerr
            << "Excluding pixels from background of source "
            << "(x = " << x << ", y = " << y << ") "
            << "bounded within aperture = " << __aperture
            << ", checking pixels with "
            << xi_min << " <= x < " << xi_max
            << ", "
            << yi_min << " <= y < " << yi_max
            << std::endl;
#endif
        for (unsigned long yi = yi_min; yi < yi_max; ++yi) {
            double y_pix_center = yi - ycent;
            for (unsigned long xi = xi_min; xi < xi_max; ++xi) {
                double x_pix_center = xi - xcent;
                if(
                    x_pix_center * x_pix_center
                    +
                    y_pix_center * y_pix_center <= ap2
                ) {
                    _bg_values(xi, yi) = Core::NaN;
                }
            }
        }
    }

} //End Background namespace.
