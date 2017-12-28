/**\file
 *
 * \brief The definitions of the functinos from CInterface.h
 *
 * \ingroup IO
 */

#define BUILDING_LIBRARY
#include "CInterface.h"

void parse_hat_mask(const char *mask_string,
                    long x_resolution,
                    long y_resolution,
                    char *mask)
{
    return IO::parse_hat_mask(mask_string, x_resolution, y_resolution, mask);
}
