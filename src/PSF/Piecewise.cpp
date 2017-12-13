/**\file
 * \brief Defines some of the Piecewise methods.
 */

#include "Piecewise.h"
#include <iostream>

namespace PSF {

    Core::vector_size_type Piecewise::cell_index(
        const std::vector<double> &grid,
        const double &x,
        Core::vector_size_type hint_lower,
        Core::vector_size_type hint_upper
    ) const
    {
        if(x == grid[hint_lower]) return hint_lower;
#ifdef DEBUG
        assert(x > grid[hint_lower]);
        assert(x <= grid[hint_upper]);
#endif
        std::vector<double>::const_iterator cell = std::lower_bound(
            grid.begin() + hint_lower,
            grid.begin() + hint_upper,
            x
        );
#ifdef DEBUG
        assert(cell != grid.begin());
#endif
        return cell - grid.begin() - 1;
    }

    void Piecewise::cell_span(const std::vector<double> &grid,
                                 double min_x,
                                 double max_x,
                                 Core::vector_size_type &low, 
                                 Core::vector_size_type &high) const
    {
        low = cell_index(grid, min_x);
        Core::vector_size_type jump = 1;
        while(jump + low<grid.size() && grid[low + jump] < max_x)
            jump *= 2;
        if(jump + low >= grid.size())
            high = cell_index(grid, max_x, low + jump / 2, grid.size() - 1);
        else high = cell_index(grid, max_x, low + jump / 2, low + jump);
    }

    std::valarray<double> Piecewise::integrate_single_column_rectangle(
        Core::vector_size_type x_index,
        Core::vector_size_type first_y_index,
        Core::vector_size_type last_y_index,
        double xmin_cell,
        double xmax_cell,
        double ymin_cell,
        double ymax_cell,
        const std::vector< std::valarray<double> >& parameter_sets
    ) const
    {
        Core::vector_size_type first_y_cell = cell_index(x_index, 
                                                         first_y_index),
                               last_y_cell = cell_index(x_index,
                                                        last_y_index);
        std::valarray<double> result(parameter_sets.size()
                                     ? parameter_sets[0].size() / 16
                                     : 1);
        if(parameter_sets.size())
            result = (
                __cells[first_y_cell]->integrate_partial_vspan(
                    xmin_cell,
                    xmax_cell,
                    ymin_cell,
                    true,
                    parameter_sets[first_y_cell]
                )
                +
                __cells[last_y_cell]->integrate_partial_vspan(
                    xmin_cell,
                    xmax_cell,
                    ymax_cell,
                    false,
                    parameter_sets[last_y_cell]
                )
            );
        else result = (
            __cells[first_y_cell]->integrate_partial_vspan(xmin_cell,
                                                           xmax_cell, 
                                                           ymin_cell,
                                                           true)
            +
            __cells[last_y_cell]->integrate_partial_vspan(xmin_cell,
                                                          xmax_cell,
                                                          ymax_cell,
                                                          false)
        );
        for(
            Core::vector_size_type y_ind = first_y_index+1;
            y_ind < last_y_index;
            ++y_ind
        ) {
            Core::vector_size_type y_cell = cell_index(x_index, y_ind);
            if(parameter_sets.size())
                result += __cells[y_cell]->integrate_vspan(
                    xmin_cell,
                    xmax_cell,
                    parameter_sets[y_cell]
                );
            else
                result += __cells[y_cell]->integrate_vspan(xmin_cell,
                                                           xmax_cell);
        }
        return result;
    }

    std::valarray<double> Piecewise::integrate_single_row_rectangle(
        Core::vector_size_type first_x_index,
        Core::vector_size_type last_x_index,
        Core::vector_size_type y_index,
        double xmin_cell,
        double xmax_cell,
        double ymin_cell,
        double ymax_cell,
        const std::vector< std::valarray<double> >& parameter_sets) const
    {
        Core::vector_size_type first_x_cell = cell_index(first_x_index,
                                                         y_index),
                               last_x_cell = cell_index(last_x_index,
                                                        y_index);
        std::valarray<double> result(parameter_sets.size() ?
                parameter_sets[0].size()/16 : 1);
        if(parameter_sets.size())
            result=__cells[first_x_cell]->integrate_partial_hspan(ymin_cell,
                    ymax_cell, xmin_cell, true, parameter_sets[first_x_cell])
                +
                __cells[last_x_cell]->integrate_partial_hspan(ymin_cell,
                        ymax_cell, xmax_cell, false,
                        parameter_sets[last_x_cell]);
        else result=__cells[first_x_cell]->integrate_partial_hspan(ymin_cell,
                    ymax_cell, xmin_cell, true)
                +
                __cells[last_x_cell]->integrate_partial_hspan(ymin_cell,
                        ymax_cell, xmax_cell, false);
        for(
            Core::vector_size_type x_ind = first_x_index + 1;
            x_ind<last_x_index;
            ++x_ind
        ) {
            Core::vector_size_type cell = cell_index(x_ind, y_index);
            if(parameter_sets.size())
                result+=__cells[cell]->integrate_hspan(ymin_cell, ymax_cell,
                        parameter_sets[cell]);
            else result+=__cells[cell]->integrate_hspan(ymin_cell, ymax_cell);
        }
        return result;
    }

    std::valarray<double> Piecewise::integrate_multi_row_column_rectangle(
        Core::vector_size_type first_x_index,
        Core::vector_size_type last_x_index,
        Core::vector_size_type first_y_index,
        Core::vector_size_type last_y_index,
        double xmin_cell,
        double xmax_cell,
        double ymin_cell,
        double ymax_cell,
        const std::vector< std::valarray<double> >& parameter_sets
    ) const
    {
        Core::vector_size_type bot_left_cell = cell_index(first_x_index,
                                                          first_y_index),
                               bot_right_cell = cell_index(last_x_index,
                                                           first_y_index),
                               top_left_cell = cell_index(first_x_index,
                                                          last_y_index),
                               top_right_cell = cell_index(last_x_index, 
                                                           last_y_index);
        std::valarray<double> result(parameter_sets.size() 
                                     ? parameter_sets[0].size()/16
                                     : 1);
        if(parameter_sets.size())
            result = (
                __cells[bot_left_cell]->integrate_rectangle(
                    xmin_cell, 
                    ymin_cell,
                    true,
                    true,
                    parameter_sets[bot_left_cell]
                )
                +
                __cells[bot_right_cell]->integrate_rectangle(
                    xmax_cell,
                    ymin_cell,
                    false,
                    true,
                    parameter_sets[bot_right_cell]
                )
                +
                __cells[top_left_cell]->integrate_rectangle(
                    xmin_cell,
                    ymax_cell,
                    true,
                    false,
                    parameter_sets[top_left_cell]
                )
                +
                __cells[top_right_cell]->integrate_rectangle(
                    xmax_cell,
                    ymax_cell,
                    false,
                    false,
                    parameter_sets[top_right_cell]
                )
            );
        else result = (
            __cells[bot_left_cell]->integrate_rectangle(xmin_cell,
                                                        ymin_cell,
                                                        true,
                                                        true)
            +
            __cells[bot_right_cell]->integrate_rectangle(xmax_cell,
                                                         ymin_cell,
                                                         false,
                                                         true)
            +
            __cells[top_left_cell]->integrate_rectangle(xmin_cell,
                                                        ymax_cell,
                                                        true,
                                                        false)
            +
            __cells[top_right_cell]->integrate_rectangle(xmax_cell,
                                                         ymax_cell,
                                                         false,
                                                         false)
        );
        for(
            Core::vector_size_type x_ind = first_x_index + 1;
            x_ind < last_x_index;
            ++x_ind
        ) {
            Core::vector_size_type bot_cell = cell_index(x_ind,
                                                         first_y_index),
                                   top_cell = cell_index(x_ind,
                                                         last_y_index);
            if(parameter_sets.size())
                result += (
                    __cells[bot_cell]->integrate_hspan(
                        ymin_cell,
                        true,
                        parameter_sets[bot_cell]
                    )
                    +
                    __cells[top_cell]->integrate_hspan(
                        ymax_cell,
                        false,
                        parameter_sets[top_cell]
                    )
                );
            else result += (
                __cells[bot_cell]->integrate_hspan(ymin_cell, true)
                +
                __cells[top_cell]->integrate_hspan(ymax_cell, false)
            );
        }
        for(
            Core::vector_size_type y_ind = first_y_index+1; 
            y_ind < last_y_index;
            ++y_ind
        ) {
            Core::vector_size_type left_cell = cell_index(first_x_index,
                                                          y_ind),
                                   right_cell = cell_index(last_x_index, 
                                                           y_ind);
            if(parameter_sets.size())
                result += (
                    __cells[left_cell]->integrate_vspan(
                        xmin_cell,
                        true,
                        parameter_sets[left_cell]
                    )
                    +
                    __cells[right_cell]->integrate_vspan(
                        xmax_cell,
                        false,
                        parameter_sets[right_cell]
                    )
                );
            else
                result += (
                    __cells[left_cell]->integrate_vspan(xmin_cell, true)
                    +
                    __cells[right_cell]->integrate_vspan(xmax_cell, false)
                );
        }
        for(
            Core::vector_size_type x_ind = first_x_index + 1; 
            x_ind < last_x_index;
            ++x_ind
        )
            for(
                Core::vector_size_type y_ind = first_y_index + 1; 
                y_ind < last_y_index;
                ++y_ind
            ) {
                Core::vector_size_type cell = cell_index(x_ind, y_ind);
                if(parameter_sets.size())
                    result += __cells[cell]->integrate(parameter_sets[cell]);
                else
                    result += __cells[cell]->integrate();
            }
        return result;
    }

    std::valarray<double> Piecewise::integrate_rectangle_parameters(
            double center_x,
            double center_y,
            double dx,
            double dy,
            const std::vector< std::valarray<double> >&parameter_sets
    ) const
    {
        double dxh = dx / 2,
               dyh = dy / 2,
               xmin = std::max(__grid_x.front(), center_x - dxh),
               ymin = std::max(__grid_y.front(), center_y - dyh),
               xmax = std::min(__grid_x.back(), center_x + dxh),
               ymax = std::min(__grid_y.back(), center_y + dyh);
        if(
            xmin > __grid_x.back()
            ||
            xmax < __grid_x.front()
            ||
            ymin > __grid_y.back()
            ||
            ymax < __grid_y.front()
        )
            return std::valarray<double>(0.0,
                                         (parameter_sets.size()
                                          ? parameter_sets[0].size()/16
                                          : 1));
        Core::vector_size_type first_x_index,
                               last_x_index,
                               first_y_index,
                               last_y_index;
        cell_span(__grid_x, xmin, xmax, first_x_index, last_x_index);
        cell_span(__grid_y, ymin, ymax, first_y_index, last_y_index);
        double xmin_cell = xmin - __grid_x[first_x_index],
               xmax_cell = xmax - __grid_x[last_x_index],
               ymin_cell = ymin - __grid_y[first_y_index],
               ymax_cell = ymax - __grid_y[last_y_index];
        if(first_x_index == last_x_index && first_y_index == last_y_index) {
            if(parameter_sets.size()) {
                Core::vector_size_type index = cell_index(first_x_index,
                                                          first_y_index);
                return __cells[index]->integrate_rectangle(
                    xmin_cell,
                    xmax_cell,
                    ymin_cell,
                    ymax_cell,
                    parameter_sets[index]
                );
            } else
                return std::valarray<double>(
                    __cells[
                        cell_index(first_x_index, first_y_index)
                    ]->integrate_rectangle(xmin_cell,
                                           xmax_cell,
                                           ymin_cell,
                                           ymax_cell),
                    1
                );
        } else if(
            first_x_index < last_x_index
            &&
            first_y_index < last_y_index
        )
            return integrate_multi_row_column_rectangle(first_x_index,
                                                        last_x_index,
                                                        first_y_index,
                                                        last_y_index,
                                                        xmin_cell,
                                                        xmax_cell,
                                                        ymin_cell,
                                                        ymax_cell,
                                                        parameter_sets);
        else if(first_x_index == last_x_index)
            return integrate_single_column_rectangle(first_x_index,
                                                     first_y_index,
                                                     last_y_index,
                                                     xmin_cell,
                                                     xmax_cell,
                                                     ymin_cell,
                                                     ymax_cell,
                                                     parameter_sets);
        else
            return integrate_single_row_rectangle(first_x_index,
                                                  last_x_index,
                                                  first_y_index,
                                                  xmin_cell,
                                                  xmax_cell,
                                                  ymin_cell,
                                                  ymax_cell,
                                                  parameter_sets);
    }

    double Piecewise::integrate_single_column_wedge(
        Core::vector_size_type x_index,
        Core::vector_size_type tip_y_index,
        Core::vector_size_type corner_y_index,
        double cell_corner_x,
        double cell_corner_y,
        double cell_tip_y,
        double radius
    ) const
    {
        bool tip_up = tip_y_index > corner_y_index;
        double circle_x = -__grid_x[x_index],
               result = __cells[
                   cell_index(x_index, tip_y_index)
               ]->integrate_hcircle_piece(cell_tip_y,
                                          cell_corner_x,
                                          radius,
                                          circle_x,
                                          -__grid_y[tip_y_index],
                                          !tip_up)
                +
                __cells[
                    cell_index(x_index,corner_y_index)
                ]->integrate_hcircle_piece(cell_corner_y,
                                           cell_corner_x,
                                           radius,
                                           circle_x,
                                           -__grid_y[corner_y_index],
                                           tip_up);
        Core::vector_size_type first_y_ind, last_y_ind;
        if(tip_up) {
            first_y_ind = corner_y_index + 1;
            last_y_ind = tip_y_index;
        }
        else {
            first_y_ind = tip_y_index + 1;
            last_y_ind = corner_y_index;
        }
        for(
            Core::vector_size_type y_ind = first_y_ind;
            y_ind < last_y_ind;
            ++y_ind
        )
            result += __cells[
                cell_index(x_index, y_ind)
            ]->integrate_hcircle_piece(cell_corner_x,
                                       radius,
                                       circle_x,
                                       -__grid_y[y_ind]);
        return result;
    }

    double Piecewise::integrate_single_row_wedge(
        Core::vector_size_type tip_x_index,
        Core::vector_size_type corner_x_index,
        Core::vector_size_type y_index,
        double cell_corner_x,
        double cell_corner_y,
        double cell_tip_x,
        double radius
    ) const
    {
        bool tip_right=tip_x_index>corner_x_index;
        double circle_y=-__grid_y[y_index],
            result=
                __cells[cell_index(tip_x_index, y_index)]->
                integrate_vcircle_piece(cell_tip_x, cell_corner_y, radius,
                        -__grid_x[tip_x_index], circle_y, !tip_right)
                +
                __cells[cell_index(corner_x_index, y_index)]->
                integrate_vcircle_piece(cell_corner_x, cell_corner_y, radius,
                        -__grid_x[corner_x_index], circle_y, tip_right);
        Core::vector_size_type first_x_ind, last_x_ind;
        if(tip_right) {first_x_ind=corner_x_index+1; last_x_ind=tip_x_index;}
        else {first_x_ind=tip_x_index+1; last_x_ind=corner_x_index;}
        for(
            Core::vector_size_type x_ind=first_x_ind;
            x_ind<last_x_ind;
            ++x_ind
        )
            result += __cells[
                cell_index(x_ind, y_index)
            ]->integrate_vcircle_piece(cell_corner_y,
                                       radius,
                                       -__grid_x[x_ind],
                                       circle_y);
        return result;
    }

    double Piecewise::integrate_inside_arc(
        Core::vector_size_type start_col,
        Core::vector_size_type end_col,
        Core::vector_size_type row,
        double radius
    ) const
    {
        int step=(start_col<end_col ? 1 : -1);
        double circle_y=-__grid_y[row];
        double result=__cells[cell_index(end_col, row)]->integrate_wedge(
                    radius, -__grid_x[end_col], circle_y);
        for(
            Core::vector_size_type col = start_col;
            col != end_col;
            col += step
        )
            result += __cells[cell_index(col, row)]->integrate_vcircle_piece(
                radius,
                -__grid_x[col],
                circle_y
            );
        return result;
    }

    double Piecewise::integrate_bottom_row(
        Core::vector_size_type corner_col,
        Core::vector_size_type corner_row,
        Core::vector_size_type crossing_col,
        Core::vector_size_type tip_col,
        double cell_crossing_x,
        double cell_tip_x,
        double cell_corner_x,
        double cell_corner_y,
        double radius,
        bool up
    ) const
    {
        int right = (tip_col > corner_col ? 1 : -1);
        double circle_y = -__grid_y[corner_row],
               result = 0;
        const PiecewiseCell* cell = __cells[cell_index(corner_col,
                                                          corner_row)];
        if(crossing_col == corner_col)
            result = cell->integrate_partial_vspan(
                std::min(cell_corner_x, cell_crossing_x),
                std::max(cell_corner_x, cell_crossing_x),
                cell_corner_y,
                up
            );
        else {
            result = cell->integrate_rectangle(cell_corner_x,
                                               cell_corner_y,
                                               right > 0,
                                               up);
            for(
                Core::vector_size_type col = corner_col + right;
                col != crossing_col;
                col += right
            )
                result += __cells[
                    cell_index(col, corner_row)
                ]->integrate_hspan(cell_corner_y, up);
            cell = __cells[cell_index(crossing_col, corner_row)];
            result += cell->integrate_rectangle(cell_crossing_x,
                                                cell_corner_y,
                                                right < 0,
                                                up);
        }
        if(crossing_col != tip_col) {
            result += cell->integrate_vcircle_piece(cell_crossing_x,
                                                    cell_corner_y,
                                                    radius,
                                                    -__grid_x[crossing_col],
                                                    circle_y,
                                                    right > 0);
            for(
                Core::vector_size_type col = crossing_col + right;
                col != tip_col;
                col += right
            )
                result += __cells[
                    cell_index(col, corner_row)
                ]->integrate_vcircle_piece(cell_corner_y,
                                           radius,
                                           -__grid_x[col],
                                           circle_y);
            result += __cells[
                cell_index(tip_col, corner_row)
            ]->integrate_vcircle_piece(cell_tip_x,
                                       cell_corner_y,
                                       radius,
                                       -__grid_x[tip_col],
                                       circle_y,
                                       right < 0);
        } else {
            result += cell->integrate_vcircle_piece(
                std::min(cell_crossing_x, cell_tip_x),
                std::max(cell_crossing_x, cell_tip_x),
                cell_corner_y,
                radius,
                -__grid_x[tip_col],
                circle_y
            );
        }
        return result;
    }

    double Piecewise::integrate_middle_row(
        Core::vector_size_type row,
        Core::vector_size_type corner_col,
        Core::vector_size_type inner_crossing_col,
        Core::vector_size_type last_col,
        double cell_corner_x,
        double inner_cell_crossing_x,
        double radius,
        int right
    ) const
    {
        const PiecewiseCell *cell = __cells[cell_index(corner_col, row)];
        double circle_x = -__grid_x[inner_crossing_col],
               circle_y = -__grid_y[row],
               result = 0;
        if(inner_crossing_col == corner_col) {
            if(inner_crossing_col == last_col)
                result += cell->integrate_hcircle_piece(cell_corner_x,
                                                        radius,
                                                        circle_x,
                                                        circle_y);
            else {
                if(right > 0)
                    result += cell->integrate_vspan(cell_corner_x,
                                                    inner_cell_crossing_x);
                else result += cell->integrate_vspan(inner_cell_crossing_x,
                                                     cell_corner_x);
            }
        } else {
            result += cell->integrate_vspan(cell_corner_x, right > 0);
            for(
                Core::vector_size_type col = corner_col + right;
                col!=inner_crossing_col;
                col += right
            ) {
                result += __cells[cell_index(col, row)]->integrate();
            }
            cell = __cells[cell_index(inner_crossing_col, row)];
            if(last_col == inner_crossing_col)
                result += cell->integrate_hcircle_piece(radius,
                                                        circle_x,
                                                        circle_y);
            else 
                result += cell->integrate_vspan(inner_cell_crossing_x,
                                                right < 0);
        }
        if(inner_crossing_col != last_col) {
            result += cell->integrate_vcircle_piece(inner_cell_crossing_x,
                                                    radius,
                                                    circle_x,
                                                    circle_y,
                                                    right > 0);
            result+=integrate_inside_arc(inner_crossing_col + right,
                                         last_col,
                                         row,
                                         radius);
        }
        return result;
    }

    double Piecewise::integrate_tip_row(
        Core::vector_size_type tip_row,
        Core::vector_size_type corner_col,
        Core::vector_size_type crossing_col,
        double cell_corner_x,
        double cell_tip_y,
        double radius,
        bool up,
        bool right
    ) const
    {
        double circle_y = -__grid_y[tip_row];
        const PiecewiseCell *cell = __cells[cell_index(corner_col, 
                                                          tip_row)];
        double result = 0;
        if(crossing_col == corner_col)
            result += cell->integrate_hcircle_piece(cell_tip_y,
                                                    cell_corner_x,
                                                    radius,
                                                    -__grid_x[corner_col], 
                                                    circle_y,
                                                    !up);
        else {
            result += cell->integrate_vcircle_piece(cell_corner_x,
                                                    radius,
                                                    -__grid_x[corner_col],
                                                    circle_y,
                                                    right);
            result += integrate_inside_arc(corner_col + (right ? 1 : -1),
                                           crossing_col,
                                           tip_row,
                                           radius);
        }
        return result;
    }

    void Piecewise::impose_boundaries(
        Core::vector_size_type &index,
        Core::vector_size_type boundary1,
        Core::vector_size_type boundary2
    ) const
    {
        if((index<boundary1 && boundary2>boundary1) ||
                (index>boundary1 && boundary2<boundary1)) index=boundary1;
        else if((index>boundary2 && boundary1<boundary2) ||
                (index<boundary2 && boundary1>boundary2)) index=boundary2;
    }

    double Piecewise::integrate_multi_row_column_wedge(
        Core::vector_size_type tip_col,
        Core::vector_size_type corner_col,
        Core::vector_size_type tip_row,
        Core::vector_size_type corner_row,
        double cell_corner_x,
        double cell_corner_y,
        double cell_tip_x,
        double cell_tip_y,
        double radius
    ) const
    {
        int right = (tip_col > corner_col ? 1 : -1),
            up = (tip_row > corner_row ? 1 : -1),
            row_correction = (up > 0 ? 1 : 0);
        Core::vector_size_type
            min_crossing_col = std::min(tip_col, corner_col),
            max_crossing_col = std::max(tip_col, corner_col);
        double r2 = std::pow(radius, 2),
               crossing_x = right * std::sqrt(
                   r2
                   -
                   std::pow(__grid_y[corner_row+row_correction], 2)
               );
        if(crossing_x < __grid_x[min_crossing_col])
            crossing_x = __grid_x[min_crossing_col];
        else if(crossing_x > __grid_x[max_crossing_col + 1])
            crossing_x = __grid_x[max_crossing_col + 1];
        Core::vector_size_type crossing_col=cell_index(__grid_x,
                                                       crossing_x,
                                                       min_crossing_col,
                                                       max_crossing_col + 1);
        impose_boundaries(crossing_col, corner_col, tip_col);
        double inner_cell_crossing_x = crossing_x - __grid_x[crossing_col];
        double result=integrate_bottom_row(corner_col,
                                           corner_row,
                                           crossing_col,
                                           tip_col,
                                           inner_cell_crossing_x,
                                           cell_tip_x,
                                           cell_corner_x,
                                           cell_corner_y,
                                           radius,
                                           up > 0);
        for(
            Core::vector_size_type row = corner_row + up;
            row != tip_row;
            row += up
        ) {
            Core::vector_size_type old_crossing_col = crossing_col;
            if(right > 0) max_crossing_col = crossing_col;
            else min_crossing_col = crossing_col;
            crossing_x = right * std::sqrt(
                r2
                -
                std::pow(__grid_y[row + row_correction], 2)
            );
            if(crossing_x < __grid_x[min_crossing_col])
                crossing_x = __grid_x[min_crossing_col];
            else if(crossing_x > __grid_x[max_crossing_col + 1])
                crossing_x = __grid_x[max_crossing_col + 1];
            crossing_col=cell_index(__grid_x,
                                    crossing_x,
                                    min_crossing_col,
                                    max_crossing_col + 1);
            impose_boundaries(crossing_col, corner_col, tip_col);
            inner_cell_crossing_x = crossing_x - __grid_x[crossing_col];
            result += integrate_middle_row(row,
                                           corner_col,
                                           crossing_col,
                                           old_crossing_col,
                                           cell_corner_x,
                                           inner_cell_crossing_x,
                                           radius,
                                           right);
        }
        result += integrate_tip_row(tip_row,
                                    corner_col,
                                    crossing_col,
                                    cell_corner_x,
                                    cell_tip_y,
                                    radius,
                                    up > 0,
                                    right > 0);
        return result;
    }

    double Piecewise::integrate_wedge(double x,
                                         double y,
                                         double radius,
                                         bool left,
                                         bool bottom) const
    {
        left = (x<0 || (x == 0 && left));
        bottom = (y<0 || (y == 0 && bottom));
        double r2 = std::pow(radius, 2),
               x2 = std::pow(x, 2),
               y2 = std::pow(y, 2);
        if(r2 < x2 + y2) return 0;
        double tip_x = (left ? -1 : 1) * std::sqrt(r2 - y2),
               tip_y = (bottom ? -1 : 1) * std::sqrt(r2 - x2),
               xmin,
               xmax,
               ymin,
               ymax;

        Core::vector_size_type tip_column = cell_index(__grid_x, tip_x),
                               tip_row = cell_index(__grid_y, tip_y),
                               corner_column = cell_index(__grid_x, x),
                               corner_row = cell_index(__grid_y, y),
                               first_column,
                               last_column,
                               first_row,
                               last_row;
        if(left) {
            xmin = tip_x;
            xmax = x;
            first_column = tip_column;
            last_column = corner_column;
        } else {
            xmin = x;
            xmax = tip_x;
            first_column = corner_column;
            last_column = tip_column;
        }
        if(bottom) {
            ymin = tip_y;
            ymax = y;
            first_row = tip_row;
            last_row = corner_row;
        } else {
            ymin = y;
            ymax = tip_y;
            first_row = corner_row;
            last_row = tip_row;
        }

        double first_x_grid = __grid_x[first_column],
               first_y_grid = __grid_y[first_row];
        cell_span(__grid_x, xmin, xmax, first_column, last_column);
        cell_span(__grid_y, ymin, ymax, first_row, last_row);
        if(left) {
            tip_column = first_column;
            corner_column = last_column;
        }
        else {
            corner_column = first_column;
            tip_column = last_column;
        }
        if(bottom) {
            tip_row = first_row;
            corner_row = last_row;
        }
        else {
            corner_row = first_row;
            tip_row = last_row;
        }

        if(first_column == last_column && first_row == last_row) 
            return __cells[
                cell_index(first_column, first_row)
            ]->integrate_hcircle_piece(ymin - first_y_grid,
                                       ymax - first_y_grid,
                                       x - first_x_grid,
                                       radius,
                                       -first_x_grid,
                                       -first_y_grid);
        else if(first_column == last_column) 
            return integrate_single_column_wedge(tip_column,
                                                 tip_row,
                                                 corner_row,
                                                 x - first_x_grid,
                                                 y - __grid_y[corner_row],
                                                 tip_y - __grid_y[tip_row], 
                                                 radius);
        else if(first_row == last_row)
            return integrate_single_row_wedge(tip_column,
                                              corner_column,
                                              tip_row,
                                              x - __grid_x[corner_column],
                                              y - first_y_grid,
                                              tip_x - __grid_x[tip_column], 
                                              radius);
        else 
            return integrate_multi_row_column_wedge(
                tip_column,
                corner_column,
                tip_row,
                corner_row,
                x - __grid_x[corner_column],
                y - __grid_y[corner_row],
                tip_x - __grid_x[tip_column],
                tip_y - __grid_y[tip_row],
                radius
            );
    }

    void Piecewise::check_inside_grid(double x, double y) const
    {
        if(x<__grid_x.front() || x>__grid_x.back()) {
            std::ostringstream msg;
            msg << "Requesting Piecewise value at x=" << x
                << ", which is outside the grid x range: "
                << __grid_x.front() << " < x < " << __grid_x.back() << ".";
            throw Error::InvalidArgument("Piecewise::check_inside_grid",
                    msg.str());
        }
        if(y<__grid_y.front() || y>__grid_y.back()) {
            std::ostringstream msg;
            msg << "Requesting Piecewise value at y=" << y
                << ", which is outside the grid y range: "
                << __grid_y.front() << " < y < " << __grid_y.back() << ".";
            throw Error::InvalidArgument("Piecewise::check_inside_grid",
                    msg.str());
        }
    }

    void Piecewise::set_cell(
        const PiecewiseCell* cell,
        const Core::vector_size_type &x_index,
        const Core::vector_size_type &y_index,
        bool clone
    )
    {
        Core::vector_size_type index=cell_index(x_index, y_index);
        __cells[index]=(clone ? cell->clone() : cell);
        __cloned_cell[index]=clone;
    }

    std::valarray<double> Piecewise::operator()(
            double x,
            double y,
            const std::vector< std::valarray<double> >& parameter_sets
    ) const
    {
        if(
                x <__grid_x.front() || x > __grid_x.back()
                ||
                y<__grid_y.front() || y>__grid_y.back()
        ) 
            return std::valarray<double>(
                    0.0,
                    (parameter_sets.size() ? parameter_sets[0].size()/16 : 1)
            );
        Core::vector_size_type x_index = cell_index(__grid_x, x),
                               y_index = cell_index(__grid_y, y);
        Core::vector_size_type index = cell_index(x_index, y_index);
        return (*__cells[index])(x - __grid_x[x_index],
                                 y - __grid_y[y_index],
                                 (parameter_sets.size()
                                  ? parameter_sets[index]
                                  : std::valarray<double>()));
    }

}
