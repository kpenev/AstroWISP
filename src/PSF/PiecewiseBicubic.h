/**\file
 *
 * \brief Defines a continuously differentiable PSF where each cell is a
 * bicubic function.
 */

#ifndef __BICUBIC_PIECEWISE_PSF_H
#define __BICUBIC_PIECEWISE_PSF_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Piecewise.h"
#include "PiecewiseBicubicCell.h"
#include <iostream>

namespace PSF {

    /**\brief A Piecewise defined as a bicubic interpolation over a grid,
     * where at each grid point the value, the x and y derivatives and the
     * cross-derivative are specified.
     */
    class LIB_PUBLIC PiecewiseBicubic : public Piecewise {
    private:
        Core::vector_size_type
            ///The number of cell in the horizontal direction.
            __x_resolution,

            ///The number of cell in the vertical direction.
            __y_resolution;

        ///The cells comprising the PSF.
        std::vector<PiecewiseBicubicCell*> __cells;
    public:
        ///\brief Create a PSF model on a specified grid with invalid cells.
        ///
        ///Use set_values to finish specifying the PSF.
        template<class RandomAccessIterator>
        PiecewiseBicubic(
            ///The first horizonatal boundary.
            RandomAccessIterator first_x,

            ///One past the last horizontal boundary.
            RandomAccessIterator last_x,

            ///The first vertical boundary.
            RandomAccessIterator first_y,

            ///One past the last vertical boundary.
            RandomAccessIterator last_y
        ) :
            Piecewise(first_x, last_x, first_y, last_y),
            __x_resolution(last_x - first_x - 1),
            __y_resolution(last_y - first_y - 1),
            __cells(__x_resolution * __y_resolution, NULL)
        {}

        ///\brief Set the PSF values and derivatives.
        template<class RandomAccessIterator>
        void set_values(
            ///The PSF value on the first grid point. Subsequent values
            ///should have the x index varying faster than the y index.
            RandomAccessIterator first_value,

            ///The first x derivative.
            RandomAccessIterator first_x_deriv,

            ///The first y derivative.
            RandomAccessIterator first_y_deriv,

            ///The first xy cross-derivative.
            RandomAccessIterator first_xy_deriv
        );

        ///Deallocate any cells that were allocated.
        ~PiecewiseBicubic()
        {
            for(size_t i = 0; i < __cells.size(); ++i)
                if(__cells[i]) delete __cells[i];
        }
    };

    ///Calculates the coefficients of the bicubic polynomial for a cell
    template<class RandomAccessIterator>
    void calc_cell_coef(
        ///The PSF value at the bottom left of the grid cell.
        RandomAccessIterator value,

        ///The x derivative at the bottom left of the grid cell.
        RandomAccessIterator x_deriv,

        ///The y derivative at the bottom left of the grid cell.
        RandomAccessIterator y_deriv,

        ///The xy cross-derivative at the bottom left of the grid cell.
        RandomAccessIterator xy_deriv,

        ///The number of PSF cells in a single row.
        Core::vector_size_type psf_x_resolution,

        ///The width of the cell.
        double width,

        ///The height of the cell.
        double height,

        ///On return, filled with the calculated coefficients.
        std::valarray<double> &cell_coef,

        ///The cell_coef array is filled from offset to offset+16.
        size_t cell_coef_offset = 0
    )
    {
        double v00 = *value,
               dx00 = *x_deriv,
               dy00 = *y_deriv,
               dxy00 = *xy_deriv,
               v10 = *(value + 1),
               dx10 = *(x_deriv + 1),
               dy10 = *(y_deriv + 1),
               dxy10 = *(xy_deriv + 1),
               v01 = *(value + psf_x_resolution + 1),
               dx01 = *(x_deriv + psf_x_resolution + 1),
               dy01 = *(y_deriv + psf_x_resolution + 1),
               dxy01 = *(xy_deriv + psf_x_resolution + 1),
               v11 = *(value + psf_x_resolution + 2),
               dx11 = *(x_deriv + psf_x_resolution + 2),
               dy11 = *(y_deriv + psf_x_resolution + 2),
               dxy11 = *(xy_deriv + psf_x_resolution + 2),
               width2 = std::pow(width, 2),
               width3 = width2 * width,
               height2 = std::pow(height, 2),
               height3 = height2 * height;

        //x0y0 - OK
        cell_coef[cell_coef_offset] = v00;

        //x1y0 - OK
        cell_coef[cell_coef_offset + 1] = dx00;

        //x2y0 - OK
        cell_coef[cell_coef_offset+2] = (3.0 * (v10 - v00) / width2
                                         -
                                         (2.0 * dx00 + dx10) / width);

        //x3y0 - OK
        cell_coef[cell_coef_offset+3] = (
            v10 - v00 - dx00*width
            -
            cell_coef[cell_coef_offset+2] * width2
        ) / width3;

        //x0y1 - OK
        cell_coef[cell_coef_offset + 4] = dy00;

        //x1y1 - OK
        cell_coef[cell_coef_offset + 5] = dxy00;

        //x2y1 - OK
        cell_coef[cell_coef_offset + 6] = (3.0 * (dy10 - dy00) / width2
                                           -
                                           (dxy10 + 2.0 * dxy00) / width);

        //x3y1 - OK
        cell_coef[cell_coef_offset + 7] = (2.0 * (dy00 - dy10) / width3
                                           +
                                           (dxy00 + dxy10) / width2);

        //x0y2 - OK
        cell_coef[cell_coef_offset + 8] = (3.0 * (v01 - v00) / height2
                                           -
                                           (2.0 * dy00 + dy01) / height);


        //x1y2 - OK
        cell_coef[cell_coef_offset + 9]=(3.0 * (dx01 - dx00) / height2
                                         -
                                         (2.0 * dxy00 + dxy01) / height);

        //x2y2 - OK
        cell_coef[cell_coef_offset+10] = (
            9.0 * (v11 - v01)
            -
            6.0 * dx01 * width
            +
            (2.0 * dxy01 + dxy11) * width * height
            -
            2.0 * cell_coef[cell_coef_offset + 6] * width2 * height
            +
            3.0 * (height * (dy01 - dy11)
                   -
                   dx11 * width
                   -
                   cell_coef[cell_coef_offset + 2]* width2)
        ) / width2 / height2;

        //x3y2 - OK
        cell_coef[cell_coef_offset + 11] = (
            (dx11 - dx10) / height2
            -
            (
                (2.0 * dxy10 + dxy11) / height
                +
                (
                    cell_coef[cell_coef_offset + 9]
                    +
                    2.0 * cell_coef[cell_coef_offset + 10] * width
                )
            ) / 3.0
        ) / width2;

        //x0y3 - OK
        cell_coef[cell_coef_offset + 12]=(
            v01
            -
            v00
            -
            dy00 * height
            -
            cell_coef[cell_coef_offset + 8] * height2
        ) / height3;

        //x1y3 - OK
        cell_coef[cell_coef_offset + 13] = (2.0 * (dx00 - dx01) / height3
                                            +
                                            (dxy00 + dxy01) / height2);

        //x2y3 - OK
        cell_coef[cell_coef_offset + 14] = (
            dy11
            -
            dy01
            -
            (
                (2.0 * dxy01 + dxy11) * width
                +
                (
                    2.0 * cell_coef[cell_coef_offset + 10] * height
                    +
                    cell_coef[cell_coef_offset + 6]
                ) * width2
            ) / 3.0
        ) / width2 / height2;

        //x3x3 - OK
        cell_coef[cell_coef_offset + 15] = (
            2.0 * (v01 - v11)
            +
            width * (dx01 + dx11)
            -
            width3 * (cell_coef[cell_coef_offset + 3]
                      +
                      cell_coef[cell_coef_offset + 7] * height
                      +
                      cell_coef[cell_coef_offset + 11] * height2)
        ) / width3 / height3;
#ifdef DEBUG
        for(size_t i = 0; i < 16; ++i)
            assert(!std::isnan(cell_coef[cell_coef_offset + i]));
#endif
    }


    template<class RandomAccessIterator>
        void PiecewiseBicubic::set_values(
            RandomAccessIterator first_value,
            RandomAccessIterator first_x_deriv,
            RandomAccessIterator first_y_deriv,
            RandomAccessIterator first_xy_deriv
        )
    {
        Core::vector_size_type cell_index = 0;
        std::valarray<double> cell_coef(16);
        for(Core::vector_size_type yind = 0; yind < __y_resolution; ++yind) {
            for(
                Core::vector_size_type xind = 0;
                xind < __x_resolution;
                ++xind
            ) {
                double cell_width = grid_column_width(xind),
                       cell_height = grid_row_height(yind);
                calc_cell_coef(first_value,
                               first_x_deriv,
                               first_y_deriv,
                               first_xy_deriv,
                               __x_resolution,
                               cell_width,
                               cell_height,
                               cell_coef);
                if(__cells[cell_index]) delete __cells[cell_index];
                __cells[cell_index] = new PiecewiseBicubicCell(
                    cell_width,
                    cell_height,
                    &(cell_coef[0]),
                    &(cell_coef[0]) + 16
                );
                ++first_value;
                ++first_x_deriv;
                ++first_y_deriv;
                ++first_xy_deriv;
                ++cell_index;
            }
            ++first_value;
            ++first_x_deriv;
            ++first_y_deriv;
            ++first_xy_deriv;
        }
        set_all_cells(__cells.begin());
    }

} //End PSF namespace.

#endif
