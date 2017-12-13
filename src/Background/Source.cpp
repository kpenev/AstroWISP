/**\file
 * 
 * \brief Defines human readable output of background sources.
 *
 * \ingroup Background
 */

#include "Source.h"

namespace Background {

    std::ostream &operator<<(std::ostream &os, const Source& source)
    {
        os << "value = " << source.value()
           << ", error = " << source.error()
           << ", pixels " << source.pixels();
        return os;
    }

}
