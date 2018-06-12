#include "parse_grid.h"

namespace IO {
    PSF::Grid parse_grid_string(const std::string &grid_string)
    {
        PSF::Grid result;

        size_t split_pos = grid_string.find_first_of(';');
        std::list<double> x_grid_list = Core::parse_real_list(
            grid_string.substr(0, split_pos),
            "grid:x",
            2,
            grid_string.size()
        );
        x_grid_list.sort();
        result.x_grid.assign(x_grid_list.begin(), x_grid_list.end());

        if(split_pos < grid_string.size() - 1) {
            std::list<double> y_grid_list = Core::parse_real_list(
                grid_string.substr(split_pos + 1, std::string::npos),
                "grid:y",
                2,
                grid_string.size()
            );
            y_grid_list.sort();
            result.y_grid.assign(y_grid_list.begin(), y_grid_list.end());
        } else result.y_grid = result.x_grid;

        return result;
    }

    std::string represent_grid(const PSF::Grid &grid)
    {
        std::ostringstream result;
        for(unsigned i = 0; i < grid.x_grid.size(); ++i) {
            result << grid.x_grid[i];
            if(i != grid.x_grid.size()-1) result << ",";
            else result << ";";
        }
        for(unsigned i = 0; i < grid.y_grid.size(); ++i) {
            result << grid.y_grid[i];
            if(i != grid.y_grid.size() - 1) result << ",";
        }
        return result.str();

    }
}
