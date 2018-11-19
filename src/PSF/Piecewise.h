/**\file
 * \brief Defines the #PSF::Piecewise PSF model.
 */

#ifndef __PIECEWISE_PSF_H
#define __PIECEWISE_PSF_H

#include "../Core/SharedLibraryExportMacros.h"
#include "PiecewiseCell.h"
#include "PSF.h"
#include "../Core/Typedefs.h"
#include <cassert>
#include <vector>
#include <valarray>
#include <algorithm>

namespace PSF {

    /**\brief A PSF model consisting of a grid of rectangles each with a
     * separate intensity function.
     */
    class LIB_PUBLIC Piecewise : public PSF {
    private:
        std::vector<double>
            ///The horizontal boundaries of the PSF grid cells.
            __grid_x,

            ///The vertical boundaries of the PSF grid cells.
            __grid_y;

        ///The cells indexed as x_ind + y_ind*Nx.
        std::vector<const PiecewiseCell*> __cells;

        ///Flags for whether the corresponding cell was cloned.
        std::vector<bool> __cloned_cell;

        ///The index of the cell with the given x and y indices.
        Core::vector_size_type cell_index(
            ///The x index of the cell to get.
            Core::vector_size_type x_ind,

            ///The y index of the cell to get.
            Core::vector_size_type y_ind
        ) const
        {return x_ind + y_ind * (__grid_x.size() - 1);}

        ///\brief The index of the cell to which the given coordinate belongs.
        ///
        ///Cell boundaries are assumed to belong to the smaller cell index,
        ///except the boundary at hint_lower which belongs to the hint_lower
        ///cell.
        Core::vector_size_type cell_index(
            ///The cell boundaries.
            const std::vector<double> &grid,

            ///The coordinate whose cell index is required.
            const double &x,

            ///An index of a cell boundary known to be less than x.
            Core::vector_size_type hint_lower,

            ///An index of a cell boundary known to be greater than x.
            Core::vector_size_type hint_upper
        ) const;

        ///\brief The index of the cell to which the given coordinate belongs.
        ///
        ///Cell boundaries are assumed to belong to the smaller cell index,
        ///except the smallest boundary which belongs to the first cell.
        Core::vector_size_type cell_index(
            ///The cell boundaries.
            const std::vector<double> &grid,

            ///The coordinate whose cell index is required.
            const double &x
        ) const
        {return cell_index(grid, x, 0, grid.size() - 1);}

        ///The index of the cell to which the given coordinate belongs.
        Core::vector_size_type cell_index(
            ///The cell boundaries.
            const std::vector<double> &grid,

            ///The coordinate whose cell index is required.
            const double &x,

            ///An index of a cell boundary known to be less than x.
            Core::vector_size_type hint_lower
        ) const
        {return cell_index(grid, x, hint_lower, grid.size());}

        ///\brief If an index is outside a range it is overwritten by the
        ///closest boundary.
        void impose_boundaries(
            ///The index to keep within bouds.
            Core::vector_size_type &index,

            ///One of the bounds to keep the index within. Need not be the lower
            ///bound.
            Core::vector_size_type boundary1,

            ///The other of the bounds to keep the index within.
            Core::vector_size_type buondary2
        ) const;

        ///Finds the range of cells that tightly cover a range of positions.
        void cell_span(
            ///The cell boundaries.
            const std::vector<double> &grid,

            ///The lower bound of the range to cover.
            double min_x,

            ///The upper bound of the range to cover.
            double max_x,

            ///On output set to the index of the first cell that covers the
            ///range.
            Core::vector_size_type &low,

            ///On output set to the index of the last cell that covers the
            ///range.
            Core::vector_size_type &high
        ) const;

        ///\brief Throws an exception if the specified location is not not 
        ///inside the grid range.
        void check_inside_grid(
            ///The x coordinate of the location to check.
            double x,

            ///The y coordinate of the location to check.
            double y
        ) const;

        ///\brief Calculates the integral over a rectangle that fits in a
        ///single column of cells.
        ///
        /// \image html PiecewisePSF_single_column_rectangle_iiidddd.png "The area being integrated along with the meaning of all function arguments."
        std::valarray<double> integrate_single_column_rectangle(
            ///The index in __grid_x of the cell containing the rectangle.
            Core::vector_size_type x_index,

            ///The index of the first y cell that overlaps with the
            ///rectangle.
            Core::vector_size_type first_y_index,

            ///The index of the last y cell that overlaps with the
            ///rectangle.
            Core::vector_size_type last_y_index,

            ///The left boundary of the rectangle in cell coordinates.
            double xmin_cell,

            ///The right boundary of the rectangle in cell coordinates.
            double xmax_cell,

            ///The bottom boundary of the rectangle in cell coordinates.
            double ymin_cell,

            ///The top boundary of the rectangle in cell coordinates.
            double ymax_cell,

            ///Optional sets of coefficients to try. The sets of coefficients
            ///to try for a cell with index i are parameter_sets[i]. Each set
            ///is 16 values long and all sets for a single pixel are
            ///concatenated.
            const std::vector< std::valarray<double> >&
            parameter_sets=std::vector< std::valarray<double> >()
        ) const;

        ///\brief Calculates the integral over a rectangle that fits in a
        ///single row of cells.
        ///
        /// \image html PiecewisePSF_single_row_rectangle_iiidddd.png "The area being integrated along with the meaning of all function arguments."
        std::valarray<double> integrate_single_row_rectangle(
            ///The index of the first x cell that overlaps with the
            ///rectangle.
            Core::vector_size_type first_x_index,

            ///The index of the last x cell that overlaps with the
            ///rectangle.
            Core::vector_size_type last_x_index,

            ///The index of the y cell that contains the rectagle.
            Core::vector_size_type y_index,

            ///The left boundary of the rectangle in cell coordinates.
            double xmin_cell,

            ///The right boundary of the rectangle in cell coordinates.
            double xmax_cell,

            ///The bottom boundary of the rectangle in cell coordinates.
            double ymin_cell,

            ///The top boundary of the rectangle in cell coordinates.
            double ymax_cell,

            ///Optional sets of coefficients to try. The sets of coefficients
            ///to try for a cell with index i are parameter_sets[i]. Each set
            ///is 16 values long and all sets for a single pixel are
            ///concatenated.
            const std::vector< std::valarray<double> >&
            parameter_sets=std::vector< std::valarray<double> >()
        ) const;

        ///\brief Calculates the integral over a rectangle spanning multiple 
        ///rows and columns.
        ///
        /// \image html PiecewisePSF_multi_row_column_rectangle_iiiidddd.png "The area being integrated along with the meaning of all function arguments."
        std::valarray<double> integrate_multi_row_column_rectangle(
            ///The index of the first x cell that overlaps with the
            ///rectangle.
            Core::vector_size_type first_x_index,

            ///The index of the last x cell that overlaps with the
            ///rectangle.
            Core::vector_size_type last_x_index,

            ///The index of the first y cell that overlaps with the
            ///rectangle.
            Core::vector_size_type first_y_index,

            ///The index of the last y cell that overlaps with the
            ///rectangle.
            Core::vector_size_type last_y_index,

            ///The left boundary of the rectangle in cell coordinates.
            double xmin_cell,

            ///The right boundary of the rectangle in cell coordinates.
            double xmax_cell,

            ///The bottom boundary of the rectangle in cell coordinates.
            double ymin_cell,

            ///The top boundary of the rectangle in cell coordinates.
            double ymax_cell,

            ///Optional sets of coefficients to try. The sets of coefficients
            ///to try for a cell with index i are parameter_sets[i]. Each set
            ///is 16 values long and all sets for a single pixel are
            ///concatenated.
            const std::vector< std::valarray<double> >&
            parameter_sets=std::vector< std::valarray<double> >()
        ) const;

        ///\brief Calculates the integral over a wedge that fits in a single 
        ///column of cells.
        ///
        ///The area being integrated matches one of the following diagrams.
        /// \image html PiecewisePSF_single_column_wedge_iiidddd_ur.png
        /// \image html PiecewisePSF_single_column_wedge_iiidddd_ul.png
        /// \image html PiecewisePSF_single_column_wedge_iiidddd_dr.png
        /// \image html PiecewisePSF_single_column_wedge_iiidddd_dl.png
        double integrate_single_column_wedge(
            ///The index in __grid_x of the cell containing the wedge.
            Core::vector_size_type x_index,

            ///The index of the cell in __grid_y containing the tip of the
            ///wedge.
            Core::vector_size_type tip_y_index,

            ///The index of the cell in __grid_y containing the corner of the
            ///wedge.
            Core::vector_size_type corner_y_index,

            ///The x coordinate of the wedge corner in cell coordinates.
            double cell_corner_x,

            ///The y coordinate of the wedge corner in cell coordinates.
            double cell_corner_y,

            ///The y coordinate of the y tip of the wedge in cell
            ///coordinates.
            double cell_tip_y,

            ///The radius of curvature of the wedge arc.
            double radius
        ) const;

        ///\brief Calculates the integral over a wedge that fits in a single row
        ///of cells.
        double integrate_single_row_wedge(
            ///The index in __grid_x of the cell containing the tip of the
            ///wedge.
            Core::vector_size_type tip_x_index,

            ///The index in __grid_x of the cell containing the corner of
            ///the wedge.
            Core::vector_size_type corner_x_index,

            ///The index of the cell in __grid_y containing the wedge.
            Core::vector_size_type y_index,

            ///The x coordinate of the wedge corner in cell coordinates.
            double cell_corner_x,

            ///The y coordinate of the wedge corner in cell coordinates.
            double cell_corner_y,

            ///The x coordinate of the y tip of the wedge in cell
            ///coordinates.
            double cell_tip_x,

            ///The radius of curvature of the wedge arc.
            double radius
        ) const;

        ///\brief Integrates over the area interior to a circle that overlaps
        ///with a row of cells.
        ///
        ///The configuration of the circle and cells must be such that the
        ///circle intersect the outside vertical wall of one end cells and
        ///one of the horizontal walls of the other end cell.
        ///
        ///In addition the arc may not go through \f$n\pi/2\f$.
        double integrate_inside_arc(
            ///The column of the end cell whose vertical wall is intersected
            ///by the arc.
            Core::vector_size_type start_col,

            ///The column of the cell whose horizontal wall is intersected by
            ///the arc.
            Core::vector_size_type end_col,

            ///The row of the cells.
            Core::vector_size_type row,

            ///The radius of the circle the arc is part of.
            double radius
        ) const;

        ///\brief Integrates the row containing the corner of a multi row and
        ///column wedge.
        double integrate_bottom_row(
            ///The column of the cell containing the corner.
            Core::vector_size_type corner_col,

            ///The row of the cells being integrated.
            Core::vector_size_type corner_row,

            ///The column of the cell where the arc enters the row of cells.
            Core::vector_size_type crossing_col,

            ///The column of the cell containing the tip furthest from the
            ///corner in the x direction.
            Core::vector_size_type tip_col,

            ///The x coordinate where the arc enters the row of cells in
            ///cell coordinates.
            double cell_crossing_x,

            ///The x coordinate of the tip of the wedge in cell coordinates.
            double cell_tip_x,

            ///The x coordinate of the corner of the wedge in cell
            ///coordinates.
            double cell_corner_x,

            ///The y coordinate of the corner of the wedge in cell
            ///coordinates.
            double cell_corner_y,

            ///The radius of the circle the wedge arc is part of.
            double radius,

            ///Wether the wedge lies above the corner.
            bool up
        ) const;

        ///\brief Integrates a row betweent the corner and the tip of a multi
        ///row and column wedge.
        double integrate_middle_row(
            ///The row being integrated.
            Core::vector_size_type row,

            ///The column of the cell containing the corner.
            Core::vector_size_type corner_col,

            ///The column where the arc first enters the row being
            ///integrated.
            Core::vector_size_type inner_crossing_col,

            ///The column furthest from the corner column which still
            ///contains a piece of the wedge.
            Core::vector_size_type last_col,

            ///The x coordinate of the wedge corner in cell coordinates.
            double cell_corner_x,

            ///The innermost x where the arc enters the row in cell
            ///coordinates.
            double inner_cell_crossing_x,

            ///The radius of the circle the wedge arc is part of.
            double radius,

            ///If the wedge extends to the right/left of the corner +/- 1.
            int right
        ) const;

        ///Integrates the row of cells containing the tip of a wedge.
        double integrate_tip_row(
            ///The row being integrated.
            Core::vector_size_type tip_row,

            ///The column of the cell containing the corner.
            Core::vector_size_type corner_col,

            ///The column where the wedge crosses the row boundary.
            Core::vector_size_type crossing_col,

            ///The x coordinate of the wedge corner in cell coordinates.
            double cell_corner_x,

            ///The y coordinate of the wedge tip in cell coordinates.
            double cell_tip_y,

            ///The radius of the circle the wedge arc is part of.
            double radius,

            ///Does the wedge extend upwards of the corner?
            bool up,

            ///Does the wedge extend rightward of the corner?
            bool right
        ) const;

        ///\brief Calculates the integral over a wedge that spans multiple
        ///rows and colunms of cells.
        double integrate_multi_row_column_wedge(
            ///The index in __grid_x of the cell containing the most
            ///horizontally displaced from the corner tip of the wedge.
            Core::vector_size_type tip_col,

            ///The index in __grid_x of the cell containing the corner of
            ///the wedge.
            Core::vector_size_type corner_col,

            ///The index of the cell in __grid_y containing the most
            ///vertically displaced from the corner tip of the wedge.
            Core::vector_size_type tip_row,

            ///The index in __grid_y of the cell containing the corner of
            ///the wedge.
            Core::vector_size_type corner_row,

            ///The x coordinate of the wedge corner in cell coordinates.
            double cell_corner_x,

            ///The y coordinate of the wedge corner in cell coordinates.
            double cell_corner_y,

            ///The x coordinate of the y tip of the wedge in cell
            ///coordinates.
            double cell_tip_x,

            ///The y coordinate of the y tip of the wedge in cell
            ///coordinates.
            double cell_tip_y,

            ///The radius of curvature of the wedge arc.
            double radius
        ) const;
    protected:
        ///\brief Equivalent to integrate_rectangle_parameters() but for the
        ///current PSF setup.
        double integrate_rectangle(
            ///See integrate_rectangle_parameters()
            double center_x,

            ///See integrate_rectangle_parameters()
            double center_y,

            ///See integrate_rectangle_parameters()
            double dx,

            ///See integrate_rectangle_parameters()
            double dy
        ) const
        {
            return integrate_rectangle_parameters(center_x,
                                                  center_y,
                                                  dx,
                                                  dy)[0];
        }

        ///Integrates the PSF over the smallest of the four wedges with the
        //following boundaries:
        /// * the line x=x
        /// * the line y=y
        /// * the circle centered at (0, 0) with a radius=radius
        ///If x is 0 the left vs right wedge is chosen according to left
        ///Same for y0 and bottom.
        double integrate_wedge(
            ///The vertical boundary of the wedge.
            double x,

            ///The horizontal boundary of the wedge
            double y,

            ///The radius of the rounded boundary of the wedge.
            double radius,

            ///If x is exactly 0, this argument determines if the wedge is
            //assumed to be on the left (true) or the right(fals) side of x=0.
            bool left = false,

            ///If y is exactly 0, this argument determines if the wedge is
            //assumed to be on the bottom (true) or the top(fals) side of y=0.
            bool bottom = false
        ) const;
    public:
        ///\brief Create a PSF model on a specified grid with invalid cells.
        ///
        ///Use set_cell or set_all_cells to specify the cells.
        template<class RandomAccessIterator>
        Piecewise(
            ///The first horizonatal boundary.
            RandomAccessIterator first_x,

            ///One past the last horizontal boundary.
            RandomAccessIterator last_x,

            ///The first vertical boundary.
            RandomAccessIterator first_y,

            ///One past the last vertical boundary.
            RandomAccessIterator last_y
        ) :
            __grid_x(first_x, last_x),
            __grid_y(first_y, last_y),
            __cells((__grid_x.size() - 1) * (__grid_y.size() - 1)),
            __cloned_cell(__cells.size(), false)
        {}

        ///\brief Set one of the cells comprising the PSF.
        ///
        ///Unless the clone argument is true, the cell must not be
        ///deallocated until the Picewise object is.
        void set_cell(
                ///The new cell to store.
                const PiecewiseCell* cell,

                ///The horizontal index of the grid cell to set.
                const Core::vector_size_type &x_index,

                ///The vertical index of the grid cell to set.
                const Core::vector_size_type &y_index,

                ///Wether to store the cell directly or to clone it. The clones
                ///are destroyed by the destructor.
                bool clone=false);

        ///\brief Set all cells.
        template<class ConstCellIterator>
            void set_all_cells(
                ///An iterator over (const) pointers to cells pointing to the
                ///first cell (the one with x index = y index = 0.
                ConstCellIterator first_cell,

                ///Whether to store clones of the cells instead of the cells
                //directly. Clones are deallocated by the destructor.
                bool clone=false)
            {
                if(clone) {
                    for(Core::vector_size_type i=0; i<__cells.size(); ++i) {
                        __cells[i]=(*first_cell)->clone();
                        ++first_cell;
                    }
                } else __cells.assign(first_cell, first_cell+__cells.size());
            }

        ///Evaluates the PSF at the given position
        double operator()(
            ///The x coolrdinate to evaluate the PSF at.
            double x,

            ///The y coolrdinate to evaluate the PSF at.
            double y
        ) const
        {return operator()(x, y, std::vector< std::valarray<double> >())[0];}

        ///\brief Evaluates the PSF at the given position assuming a list of cell
        ///parameters.
        std::valarray<double> operator()(
                ///The x coordinate where to evaluate the PSF.
                double x,

                ///The y coordinate where to evaluate the PSF.
                double y,

                ///The sets of coefficients to use. The sets of coefficients
                ///to use for a cell with index i are parameter_sets[i]. Each set
                ///is 16 values long and all sets for a single pixel are
                ///concatenated. If this argument is empty, the current PSF setup
                ///is used.
                const std::vector< std::valarray<double> >& parameter_sets
        ) const;

        ///Returns the width of a cell column.
        double grid_column_width(
            ///The column to return the width of.
            Core::vector_size_type column
        ) const
        {return __grid_x[column + 1] - __grid_x[column];}

        ///Returns the height of a cell row.
        double grid_row_height(
            ///The row to return the height of.
            Core::vector_size_type row
        ) const
        {return __grid_y[row + 1] - __grid_y[row];}

        ///The left boundary of the PSF grid.
        double min_x() const {return __grid_x.front();}

        ///The right boundary of the PSF grid.
        double max_x() const {return __grid_x.back();}

        ///The bottom boundary of the PSF grid.
        double min_y() const {return __grid_y.front();}

        ///The top boundary of the PSF grid.
        double max_y() const {return __grid_y.back();}

        ///\brief Calculates the integral of the PSF over a rectangle.
        std::valarray<double> integrate_rectangle_parameters(
            ///The x coordinate of the center of the rectangle
            double center_x,

            ///The x coordinate of the center of the rectangle
            double center_y,

            ///The full width of the rectangle.
            double dx,

            ///The full height of the rectangle.
            double dy,

            ///Optional sets of coefficients to use. The sets of coefficients
            ///to use for a cell with index i are parameter_sets[i]. Each set
            ///is 16 values long and all sets for a single pixel are
            ///concatenated. If this argument is omitted or empty, the
            ///current PSF setup is used.
            const std::vector< std::valarray<double> >&
            parameter_sets=std::vector< std::valarray<double> >()
        ) const;
    };

} //End PSF namespace.

#endif
