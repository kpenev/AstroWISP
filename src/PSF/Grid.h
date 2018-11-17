/**\file
 *
 * \brief Some useful typedef statements.
 */

#ifndef __PSFGRID_H
#define __PSFGRID_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/ParseCSV.h"
#include <string>
#include <vector>
#include <list>
#include <sstream>

namespace PSF {

    ///The cell grid over which the PSF/PRF is defined.
    class LIB_PUBLIC Grid {
    public:
        ///Default constructor (grid with no cells).
        Grid() {}

        std::vector<double>
            ///The x boundaries between grid cells.
            x_grid,

            ///The y boundaries between grid cells.
            y_grid;
    };

}

#endif
