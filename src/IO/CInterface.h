/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * IO library.
 *
 * \ingroup IO
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "parse_hat_mask.h"

extern "C" {

    ///C-binding alias for IO::parse_hat_mask.
    LIB_PUBLIC char *parse_hat_mask(const char *mask_string,
                                    long x_resolution,
                                    long y_resolution);

} //End Extern "C".
