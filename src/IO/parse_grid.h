/**\file
 *
 * \brief Define functions for converting between string and PSF::Grid.
 *
 * \ingroup IO
 */


#include "../PSF/Grid.h"

#include <string>

namespace IO {
    ///Return the grid represented in the given string.
    PSF::Grid LIB_PUBLIC parse_grid_string(const std::string &grid_string);

    ///Return the string representing the given grid.
    std::string LIB_PUBLIC represent_grid(const PSF::Grid &grid);
}//End IO namespace.
