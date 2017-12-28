/**\file
 *
 * \brief Declares a function for parsing HAT-style maskes from image header.
 *
 * \ingroup IO
 */

#include "../Core/SharedLibraryExportMacros.h"

namespace IO {
    ///\brief Fill a pre-allocated array with resolution of (x_resolution x
    ///y_resolution) with the parsed mask.
    LIB_PUBLIC void parse_hat_mask(
        ///The concatenated MASKINFO entries from the header.
        const char *mask_string,

        ///The horizontal resolution of the image.
        long x_resolution,

        ///The vertical resolution of the image.
        long y_resolution,

        ///The array to write the mask to. Only the masked entries are updated
        ///using logical or, thus calling this function with different 
        ///mask_string values but the same mask array combines all masks.
        char *mask
    );

} //End IO namespace.
