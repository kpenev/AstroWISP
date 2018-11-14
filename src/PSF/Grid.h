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

    class LIB_PUBLIC Grid {
    public:
        ///Default constructor (grid with no cells).
        Grid() {}

        std::vector<double> x_grid, y_grid;
    };

}

#endif
