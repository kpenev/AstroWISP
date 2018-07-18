/**\file
 *
 * \brief Implementation for some methods of Manual.
 *
 * \ingroup Background
 */

#include "Manual.h"

namespace Background {
    Manual::Manual(double *value,
                   double *error,
                   unsigned *num_pixels,
                   unsigned num_sources)
    {
        __sources.reserve(num_sources);
        for(unsigned src_index = 0; src_index < num_sources; ++src_index)
            __sources.push_back(Source(value[src_index],
                                       error[src_index],
                                       num_pixels[src_index]));
    }

}//End Background namespace.
