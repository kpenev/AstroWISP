/**\file
 *
 * \brief Declares a smoothly varying over an image piecewise bicubic PSF map
 * class.
 */

#ifndef __PIECEWISE_BICUBIC_PSF_MAP_H
#define __PIECEWISE_BICUBIC_PSF_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../IO/parse_grid.h"
#include "PiecewiseBicubic.h"
#include "Map.h"
#include "Eigen/Dense"

#include <iostream>
#include "../Core/NaN.h"

namespace PSF {

    ///Simplify statements requiring this type.
    typedef Eigen::Map< const Eigen::Matrix<double,
                                            Eigen::Dynamic,
                                            Eigen::Dynamic,
                                            Eigen::RowMajor> > RowMajorMap;

    ///Smoothly varying over an image piecewise bicubic PSF.
    class LIB_PUBLIC PiecewiseBicubicMap : public Map {
    private:
        std::vector<double>
            ///The horizontal boundaries between PSF grid cells.
            __x_grid,

            ///The vertical boundaries between PSF grid cells.
            __y_grid;

        ///A direct pointer to the coefficients defining the map.
        const double *__coefficients;

        Eigen::MatrixXd
            ///\brief Apply to a vector of polynomial terms to get the values
            ///of the PSF at the grid intersections.
            __value_expansion,

            ///\brief Apply to a vector of polynomial terms to get the x
            ///derivatives of the PSF at the grid intersections.
            __x_deriv_expansion,

            ///\brief Apply to a vector of polynomial terms to get the y
            ///derivatives of the PSF at the grid intersections.
            __y_deriv_expansion,

            ///\brief Apply to a vector of polynomial terms to get the x-y
            ///cross derivatives of the PSF at the grid intersections.
            __xy_deriv_expansion;

        ///Was the x grid expanded to the left (1: yes, 0: no)?
        short __expanded_left,

             ///Was the x grid expanded to the right  (1: yes, 0: no)?
             __expanded_right,

             ///Was the y grid expanded upward  (1: yes, 0: no)?
             __expanded_up,

             ///Was the y grid expanded downward  (1: yes, 0: no)?
             __expanded_down;

        ///Fills the coefficients for a vector of values.
        template<class VECTOR_TYPE>
        void fill_coefficients(
                ///The values of the coefficients as they are returned by PSF
                ///fitting.
                const VECTOR_TYPE& coefficient_vector);
    public:
        ///Construct a map from the given data tree.
        PiecewiseBicubicMap(
            ///The I/O data tree to construct the map from. It should contain
            ///at least the quantities identified by
            ///PiecewiseBicubicMap::required_data_tree_quantities().
            const IO::H5IODataTree &data,

            ///The PSF grid should be expanded to cover at least a circle of
            ///this radius around each source.
            double min_psf_span = 0
        );

        ///\brief Create a PSF map from a given polynomial expansion of the
        ///PSF parameters.
        PiecewiseBicubicMap(
            ///The coefficients of the polynomial expansion of the PSF
            ///parameters. Should consist of 4 parts of equal length, giving
            ///the expansion of the values, x derivatives, y derivatives and
            ///xy cross-derivatives at each grid point. In each part, the
            ///polynomial expansion coefficients for a single grid
            ///intersection are sequential in memory, with each row of grid
            ///intersections occupying a contiguous block.
            const Eigen::VectorXd &expansion_coef,

            ///The x coordinates of the PSF cell boundaries.
            const std::vector<double> &x_grid,

            ///The y coordinates of the PSF cell boundaries.
            const std::vector<double> &y_grid
        );

        ///The horizontal boundaries between PSF grid cells.
        const std::vector<double> &x_grid() const {return __x_grid;}

        ///The vertical boundaries between PSF grid cells.
        const std::vector<double> &y_grid() const {return __y_grid;}

        ///The underlying PSF fit coefficients (as given by PSF fitting)
        const double *coefficients() const {return __coefficients;}

        ///\brief Evaluate the PSF map for a particular set of values of the
        ///terms and return the PSF.
        PiecewiseBicubic *get_psf(
            ///The values of the terms on which the PSF map depends.
            const Eigen::VectorXd &terms,

            ///Background to add to the PSF (assumed constant) normalized in
            ///the same way as the backgroundless PSFs produced by the map.
            double background = 0
        ) const;

        ///A reference to a dynamically allocated PSF.
        PiecewiseBicubic *operator()(
            ///The values of the terms on which the PSF map depends.
            const Eigen::VectorXd &terms,

            ///Background to add to the PSF (assumed constant) normalized in
            ///the same way as the backgroundless PSFs produced by the map.
            double background = 0
        ) const
        {return get_psf(terms, background);}

        ///\brief All quantities needed to construct the PSF map from an I/O
        ///data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    }; //End PiecewiseBicubicMap class.

} //End PSF namespace.
#endif
