/**\file
 *
 * \brief The definitions of the functinos from CInterface.h
 *
 * \ingroup IO
 */

#define BUILDING_LIBRARY
#include "CInterface.h"

const char MASK_OK = Core::MASK_OK;
const char MASK_CLEAR = Core::MASK_CLEAR;
const char MASK_FAULT = Core::MASK_FAULT;
const char MASK_HOT = Core::MASK_HOT;
const char MASK_COSMIC = Core::MASK_COSMIC;
const char MASK_OUTER = Core::MASK_OUTER;
const char MASK_OVERSATURATED = Core::MASK_OVERSATURATED;
const char MASK_LEAKED = Core::MASK_LEAKED;
const char MASK_SATURATED = Core::MASK_SATURATED;
const char MASK_INTERPOLATED = Core::MASK_INTERPOLATED;
const char MASK_BAD = Core::MASK_BAD;
const char MASK_ALL = Core::MASK_ALL;
const char MASK_NAN = Core::MASK_NAN;

void parse_hat_mask(const char *mask_string,
                    long x_resolution,
                    long y_resolution,
                    char *mask)
{
    return IO::parse_hat_mask(mask_string, x_resolution, y_resolution, mask);
}
