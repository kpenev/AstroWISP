/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * IO library.
 *
 * \ingroup IO
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Image.h"
#include "parse_hat_mask.h"

extern "C" {

    ///See Core::MASK_OK
    LIB_PUBLIC extern const char MASK_OK;

    ///See Core::MASK_CLEAR
    LIB_PUBLIC extern const char MASK_CLEAR;


    ///See Core::MASK_FAULT
    LIB_PUBLIC extern const char MASK_FAULT;

    ///See Core::MASK_HOT
    LIB_PUBLIC extern const char MASK_HOT;

    ///See Core::MASK_COSMIC
    LIB_PUBLIC extern const char MASK_COSMIC;

    ///See Core::MASK_OUTER
    LIB_PUBLIC extern const char MASK_OUTER;

    ///See Core::MASK_OVERSATURATED
    LIB_PUBLIC extern const char MASK_OVERSATURATED;

    ///See Core::MASK_LEAKED
    LIB_PUBLIC extern const char MASK_LEAKED;

    ///See Core::MASK_SATURATED
    LIB_PUBLIC extern const char MASK_SATURATED;

    ///See Core::MASK_INTERPOLATED
    LIB_PUBLIC extern const char MASK_INTERPOLATED;

    ///See Core::MASK_BAD
    LIB_PUBLIC extern const char MASK_BAD;

    ///See Core::MASK_ALL
    LIB_PUBLIC extern const char MASK_ALL;

    ///See Core::MASK_NAN
    LIB_PUBLIC extern const char MASK_NAN;

    ///C-binding alias for IO::parse_hat_mask.
    LIB_PUBLIC void parse_hat_mask(const char *mask_string,
                                   long x_resolution,
                                   long y_resolution,
                                   char *mask);

} //End Extern "C".
