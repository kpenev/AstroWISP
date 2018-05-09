/**\file
 *
 * \brief The implementation of the methods of the MeasureAnnulus class.
 *
 * \ingroup Background
 */

#include "MeasureAnnulus.h"

namespace Background {

    Source MeasureAnnulus::operator()(double x, double y) const
    {
        const double xoff = x - 0.5;
        const double yoff = y - 0.5;
        double rin2 = __inner_radius * __inner_radius, 
               rout2 = __outer_radius * __outer_radius;
        unsigned long
            xi_min = std::max(0.0, std::ceil(xoff - __outer_radius)),
            xi_max = std::min(static_cast<double>(_bg_values.x_resolution()), 
                            std::ceil( xoff + __outer_radius )),
            yi_min = std::max(0.0, std::floor( yoff - __outer_radius )),
            yi_max = std::min(static_cast<double>(_bg_values.y_resolution()), 
                            std::ceil( yoff + __outer_radius ));
        size_t max_bg_values = (
            std::pow(std::floor(2.0 * __outer_radius + 1.0), 2)
            -
            std::pow(std::floor(std::sqrt(2.0) * __inner_radius - 1.0), 2)
        );
        std::valarray<double> surviving_bg_values(max_bg_values);

        double* first_bg_value = &(surviving_bg_values[0]);
        double* last_bg_value = first_bg_value;

        for ( unsigned long yi = yi_min; yi < yi_max; yi++ ) {
            double y_pix_center = yi - yoff;
            double y_pix_center2 = y_pix_center * y_pix_center;

            for ( unsigned long xi = xi_min; xi < xi_max; xi++ ) {
                double x_pix_center = xi - xoff;
                double pic_dist = (
                    x_pix_center * x_pix_center + y_pix_center2
                );

                if ( pic_dist > rin2 && pic_dist <= rout2 ) {
                    double bg_val = _bg_values( xi, yi );
                    if( !std::isnan( bg_val )) {
                        *(last_bg_value++) = bg_val;
                    }
                }
            }
        }

        size_t num_bg_values = last_bg_value - first_bg_value;
        if(num_bg_values < 3) return Core::NaN;
        std::sort(first_bg_value, last_bg_value);
        double *median_bg_value = first_bg_value + num_bg_values / 2, 
               *bg_range = first_bg_value;
        size_t range_len = static_cast<size_t>(
            __error_confidence * num_bg_values + 0.5
        );
        double medbg = (num_bg_values % 2
                        ? *median_bg_value
                        : 0.5 * (*median_bg_value + *(median_bg_value - 1)));
        while(
            bg_range != last_bg_value - range_len
            &&
            medbg - bg_range[0] > bg_range[range_len] - medbg
        )
            ++bg_range;
        double min_error = std::max(medbg - bg_range[0],
                                    bg_range[range_len - 1] - medbg),
               max_error;
        if(bg_range > first_bg_value && bg_range + range_len < last_bg_value)
            max_error = std::min(medbg - bg_range[-1],
                                 bg_range[range_len] - medbg);
        else if(bg_range > first_bg_value)
            max_error = medbg - bg_range[-1];
        else if(bg_range + range_len < last_bg_value)
            max_error = bg_range[range_len] - medbg;
        else
            max_error = min_error;

        Source result(
            medbg,
            (min_error + max_error) / 2.0
            *
            std::sqrt(M_PI / (2.0 * (num_bg_values - 1))),
            num_bg_values
        );
#ifdef VERBOSE_DEBUG
        std::cerr << "BG(" << x << ", " << y << "): " << result << std::endl;
#endif
        return result;
    } //End MeasureAnnulus::operator() definition.
}
