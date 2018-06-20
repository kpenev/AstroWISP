/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup Core
 */

#include "CInterface.h"
#include "Image.h"
#include "SubPixelMap.h"

CoreImage *create_core_image(unsigned long x_resolution,
                             unsigned long y_resolution,
                             double *values,
                             double *errors,
                             char *mask,
                             bool wrap)
{
    Core::Image<double> *result;
    if(wrap) {
        result = new Core::Image<double>();
        result->wrap(values,
                     mask,
                     x_resolution,
                     y_resolution,
                     errors);
    } else {
        result = new Core::Image<double>(values,
                                         mask,
                                         x_resolution,
                                         y_resolution,
                                         errors);
    }
    return reinterpret_cast<CoreImage*>(result);
}

void destroy_core_image(CoreImage *image)
{
    delete reinterpret_cast<Core::Image<double>*>(image);
}

CoreSubPixelMap *create_core_subpixel_map(unsigned long x_resolution,
                                          unsigned long y_resolution,
                                          double *sensitivities)
{
    Core::SubPixelMap *result = new Core::SubPixelMap(x_resolution,
                                                      y_resolution,
                                                      "c_interface");
    for(unsigned y = 0; y < y_resolution; ++y)
        for(unsigned x = 0; x < x_resolution; ++x) {
            (*result)(x, y) = *sensitivities;
            ++sensitivities;
        }
    return reinterpret_cast<CoreSubPixelMap*>(result);
}

void destroy_core_subpixel_map(CoreSubPixelMap *map)
{
    delete reinterpret_cast<Core::SubPixelMap*>(map);
}
