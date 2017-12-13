/**\file
 *
 * \brief Some useful typedef statements.
 */

#ifndef __PSFGRID_H
#define __PSFGRID_H

#include "../Core/ParseCSV.h"
#include <string>
#include <vector>
#include <list>
#include <sstream>

namespace PSF {

    class Grid {
    public:
        ///Default constructor (grid with no cells).
        Grid() {}

        ///Create a grid according to a string
        Grid(
            ///The string defining the grid. The format consists of comma
            ///separated list of x grid boundaries, optionally followed by
            ///';' and a comma separated list of y grid boundaries. If the y
            ///grid boundaries part is omitted, the x_grid boundaries are
            ///used.
            const std::string &grid_string
        );

        std::vector<double> x_grid, y_grid;

        operator std::string() const;
    };

}

#endif
