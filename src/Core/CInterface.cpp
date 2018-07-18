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
#ifdef VERBOSE_DEBUG
        std::cerr << "Wrapped (" << x_resolution << "x" << y_resolution
                  << ") image around: "
                  << "\tvalues = " << values
                  << "\terrors = " << errors
                  << "\tmask = " << (void*)mask
                  << std::endl;
#endif
    } else {
        result = new Core::Image<double>(values,
                                         mask,
                                         x_resolution,
                                         y_resolution,
                                         errors);
#ifdef VERBOSE_DEBUG
        std::cerr << "Created new image" << std::endl;
#endif
    }
#ifdef VERBOSE_DEBUG
        std::cerr << "New image at " << result << std::endl;
#endif
    return reinterpret_cast<CoreImage*>(result);
}

void destroy_core_image(CoreImage *image)
{
#ifdef VERBOSE_DEBUG
        std::cerr << "Destroying image at " << image << std::endl;
#endif
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
