#ifndef __PIECEWISE_BICUBIC_PSF_SMOOTHING_H
#define __PIECEWISE_BICUBIC_PSF_SMOOTHING_H

#include "../Core/SharedLibraryExportMacros.h"
#include "LinearSource.h"
#include "Eigen/Dense"
#include <vector>
#include <valarray>
#include <iostream>

namespace FitPSF {

    ///A class that applies a smoothing penalty during PSF/PRF fitting.
    class LIB_LOCAL PiecewiseBicubicPSFSmoothing {
    private:

        Eigen::MatrixXd
            ///The average of \f$ {\kappa^s}^T \kappa^s \f$ over all sources.
            __lhs_spatial_dependence,

            ///See lhs_correction().
            __lhs_matrix;


        Eigen::VectorXd
            ///The average of the PSF terms for all sources.
            __rhs_spatial_dependence,

            ///See rhs_correction().
            __rhs_vector;

        ///Is this object ready for smoothing?
        bool __prepared;

        ///\brief Fill the matrix which converts cell parameters to the
        ///smoothing penalty.
        ///
        ///This is the \f$ \mathbf{A}_{cell} \f$ in the documentation.
        void fill_cell_penalty_matrix(
            ///The width of the cell.
            double width,

            ///The hegiht of the cell.
            double height,

            ///The matrix to fill.
            Eigen::MatrixXd& penalty_matrix
        );

        ///\brief Fills the (2,2), (2,3), (3,2), (3,3) coefficients in the
        ///bicubic function over a cell.
        inline void fill_second_derivative_coefficients(
            ///A pointer to the first coefficient defining the bicubic over
            ///the cell.
            const double *first_coef,

            ///The vector to fill.
            Eigen::Vector4d &result
        );

        ///Fill the __rhs_spatial_depedence vector.
        void fill_spatial_dependence(
            ///See same name argument to constructor.
            const LinearSourceList &fit_sources
        );

        ///\brief Updates RHS vector with the smoothing terms for single
        ///cell/basis vector.
        void add_smoothing_rhs(
            ///The basis parameter sets returned by fill_basis_parameter_sets.
            const std::vector< std::valarray<double> > &basis_parameter_sets,

            ///The penalty matrix returned by fill_cell_penalty_matrix.
            const Eigen::MatrixXd &penalty_matrix,

            ///The index of the basis vector to add (zero should correspond to
            //the unity integral constraint, so this should be >=1).
            unsigned basis_index,

            ///The index of the cell (assuming the order in
            ///basis_parameter_sets).
            unsigned cell_index,

            ///The RHS vector to add to.
            Eigen::VectorXd &rhs
        );

        ///\brief Updates LHS matrix with the smoothing terms for single
        ///cell/basis vector.
        void add_smoothing_matrix(
            ///The basis parameter sets returned by fill_basis_parameter_sets.
            const std::vector< std::valarray<double> > &basis_parameter_sets,

            ///The penalty matrix returned by fill_cell_penalty_matrix.
            const Eigen::MatrixXd &penalty_matrix,

            ///The index of the basis vector to add (zero should correspond to
            //the unity integral constraint, so this should be >=1).
            unsigned basis_index,

            ///The index of the cell (assuming the order in
            ///basis_parameter_sets).
            unsigned cell_index,

            ///The RHS vector to add to.
            Eigen::MatrixXd &matrix
        );

    public:
        ///Default constructor, must call prepare_smoothing before use.
        PiecewiseBicubicPSFSmoothing() : __prepared(false) {}

        ///Get ready for imposing smoothing penalty to bicubic PSF fits.
        PiecewiseBicubicPSFSmoothing(
            ///The degree of smoothing to use.
            double smoothing_penalty,

            ///The sources which participate in the PSF shape fitting.
            const LinearSourceList &fit_sources,

            ///The list of vertical grid cell boundaries of the PSF.
            const std::vector<double> &x_grid,

            ///The list of horizontal grid cell boundaries of the PSF.
            const std::vector<double> &y_grid,

            ///The parameter sets in the basis in which PSF fitting is
            ///performed.
            const std::vector< std::valarray<double> > &basis_parameter_sets
        )
        {prepare_smoothing(smoothing_penalty,
                           fit_sources,
                           x_grid, y_grid,
                           basis_parameter_sets);}

        ///See non-default constructor.
        void prepare_smoothing(
            double smoothing_penalty,
            const LinearSourceList &fit_sources,
            const std::vector<double> &x_grid,
            const std::vector<double> &y_grid,
            const std::vector< std::valarray<double> > &basis_parameter_sets
        );

        ///The vector to add to the RHS vector to implement smoothing.
        const Eigen::VectorXd &rhs_correction() const {return __rhs_vector;}

        ///The matrix to add to the LHS matrix to implement smoothing.
        const Eigen::MatrixXd &lhs_correction() const {return __lhs_matrix;}

        ///Is smoothing ready.
        operator bool() const {return __prepared;}
    }; //End PiecewiseBicubicPSFSmoothing class.

}//End FitPSF namespace.

#endif
