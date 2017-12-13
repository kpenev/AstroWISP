#include "Grid.h"

namespace PSF {

    Grid::Grid(const std::string &grid_string)
    {
        size_t split_pos = grid_string.find_first_of(';');
        std::list<double> x_grid_list = Core::parse_real_list(
            grid_string.substr(0, split_pos),
            "grid:x",
            2,
            grid_string.size()
        );
        x_grid_list.sort();
        x_grid.assign(x_grid_list.begin(), x_grid_list.end());
        if(split_pos < grid_string.size() - 1) {
            std::list<double> y_grid_list = Core::parse_real_list(
                grid_string.substr(split_pos + 1, std::string::npos),
                "grid:y",
                2,
                grid_string.size()
            );
            y_grid_list.sort();
            y_grid.assign(y_grid_list.begin(), y_grid_list.end());
        } else y_grid = x_grid;
    }

    Grid::operator std::string() const
    {
        std::ostringstream result;
        for(unsigned i = 0; i < x_grid.size(); ++i) {
            result << x_grid[i];
            if(i != x_grid.size()-1) result << ",";
            else result << ";";
        }
        for(unsigned i = 0; i < y_grid.size(); ++i) {
            result << y_grid[i];
            if(i != y_grid.size() - 1) result << ",";
        }
        return result.str();
    }

}
