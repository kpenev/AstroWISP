/**\file
 *
 * \brief Declares a function for parsing HAT-style maskes from image header.
 *
 * \ingroup IO
 */

#include "../Core/SharedLibraryExportMacros.h"

namespace IO {
    ///\brief Return a newly allocated array with resolution of (x_resolution x
    ///y_resolution) containing the parsed mask.
    LIB_PUBLIC char *parse_hat_mask(const char *mask_string,
                                    long x_resolution,
                                    long y_resolution);

} //End IO namespace.
