#include "PiecewiseBicubicPSFSmoothing.h"

namespace FitPSF {

    void PiecewiseBicubicPSFSmoothing::fill_cell_penalty_matrix(
            double width,
            double height,
            Eigen::MatrixXd& penalty_matrix
    )
    {
        double width2 = width * width,
               height2 = height * height;

        penalty_matrix << 16.0,
                          24.0 * height,
                          24.0 * width,
                          36.0 * width * height,

                          24.0 * height,
                          48.0 * height2,
                          36.0 * width * height,
                          72.0 * width * height2,

                          24.0 * width,
                          36.0 * width * height,
                          48.0 * width2,
                          72.0 * width2 * height,

                          36.0 * width * height,
                          72.0 * width * height2,
                          72.0 * width2 * height,
                          144.0 * width2 * height2;

        penalty_matrix *= width * height;
    }

    inline void 
    PiecewiseBicubicPSFSmoothing::fill_second_derivative_coefficients(
            const double *first_coef,
            Eigen::Vector4d &result
    )
    {
        result[0] = first_coef[10];
        result[1] = first_coef[14];
        result[2] = first_coef[11];
        result[3] = first_coef[15];
    }

    void PiecewiseBicubicPSFSmoothing::fill_spatial_dependence(
        const LinearSourceList &fit_sources
    )
    {
        LinearSourceList::const_iterator src_i = fit_sources.begin();

        __rhs_spatial_dependence = (*src_i)->expansion_terms();
        __lhs_spatial_dependence = (*src_i)->expansion_terms()
                                   *
                                   (*src_i)->expansion_terms().transpose();

        while(src_i != fit_sources.end()) {
            __rhs_spatial_dependence += (*src_i)->expansion_terms();

            __lhs_spatial_dependence += (
                (*src_i)->expansion_terms()
                *
                (*src_i)->expansion_terms().transpose()
            );

            ++src_i;
        }

        __rhs_spatial_dependence /= fit_sources.size();
        __lhs_spatial_dependence /= fit_sources.size();
    }

    void PiecewiseBicubicPSFSmoothing::add_smoothing_rhs(
            const std::vector< std::valarray<double> > &basis_parameter_sets,
            const Eigen::MatrixXd &penalty_matrix,
            unsigned basis_index,
            unsigned cell_index,
            Eigen::VectorXd &rhs
    )
    {
        assert(basis_index >= 1);
        Eigen::Vector4d basis_coef, total_integral_coef;
        fill_second_derivative_coefficients(
                &(basis_parameter_sets[cell_index][16*basis_index]),
                basis_coef
        );
        fill_second_derivative_coefficients(
                &(basis_parameter_sets[cell_index][0]),
                total_integral_coef
        );

        rhs.segment((basis_index - 1) * __rhs_spatial_dependence.size(),
                    __rhs_spatial_dependence.size()) -=
            (basis_coef.transpose() * penalty_matrix * total_integral_coef)[0]
            *
            __rhs_spatial_dependence;
    }

    void PiecewiseBicubicPSFSmoothing::add_smoothing_matrix(
            const std::vector< std::valarray<double> > &basis_parameter_sets,
            const Eigen::MatrixXd &penalty_matrix,
            unsigned row_basis_index,
            unsigned cell_index,
            Eigen::MatrixXd &matrix
    )
    {
        assert(row_basis_index >= 1);
        unsigned num_spatial_coef = __lhs_spatial_dependence.rows(),
                 num_basis_vectors = basis_parameter_sets[0].size() / 16;
        assert(__lhs_spatial_dependence.cols() == num_spatial_coef);
        Eigen::Vector4d row_basis_coef;
        fill_second_derivative_coefficients(
            &(basis_parameter_sets[cell_index][16 * row_basis_index]),
            row_basis_coef
        );
        for(
            unsigned col_basis_index = 1;
            col_basis_index < num_basis_vectors;
            ++col_basis_index
        ) {
            Eigen::Vector4d col_basis_coef;
            fill_second_derivative_coefficients(
                &(basis_parameter_sets[cell_index][16 * col_basis_index]),
                col_basis_coef
            );
            matrix.block(
                (row_basis_index - 1) * num_spatial_coef,
                (col_basis_index - 1) * num_spatial_coef,
                num_spatial_coef,
                num_spatial_coef
            ) += (
                (
                    row_basis_coef.transpose()
                    *
                    penalty_matrix
                    *
                    col_basis_coef
                )[0]
                *
                __lhs_spatial_dependence
            );
        }
    }

    void PiecewiseBicubicPSFSmoothing::prepare_smoothing(
        double smoothing_penalty,
        const LinearSourceList &fit_sources,
        const std::vector<double> &x_grid,
        const std::vector<double> &y_grid,
        const std::vector< std::valarray<double> > &basis_parameter_sets
    )
    {
        fill_spatial_dependence(fit_sources);

        size_t cell_ind = 0,
               num_basis_vectors = basis_parameter_sets[0].size() / 16.0,
               num_expansion_terms = 
                   fit_sources.front()->expansion_terms().size(),
               num_parameters = (num_basis_vectors - 1) * num_expansion_terms;

        __rhs_vector = Eigen::VectorXd::Zero(num_parameters);
        __lhs_matrix = Eigen::MatrixXd::Zero(num_parameters, num_parameters);
        for( size_t y_ind=0; y_ind<y_grid.size()-1; ++y_ind ) {
            double height=y_grid[y_ind+1]-y_grid[y_ind];
            for( size_t x_ind=0; x_ind<x_grid.size()-1; ++x_ind) {
                double width=x_grid[x_ind+1]-x_grid[x_ind];
                Eigen::MatrixXd penalty_matrix(4, 4);
                fill_cell_penalty_matrix(width, height, penalty_matrix);
                for( 
                    unsigned basis_ind=1; 
                    basis_ind<num_basis_vectors;
                    ++basis_ind
                ) {
                    add_smoothing_rhs(basis_parameter_sets,
                                      penalty_matrix,
                                      basis_ind,
                                      cell_ind,
                                      __rhs_vector);
                    add_smoothing_matrix(basis_parameter_sets,
                                         penalty_matrix,
                                         basis_ind,
                                         cell_ind,
                                         __lhs_matrix);
                }
                ++cell_ind;
            }
        }

        __rhs_vector *= std::pow(10.0, smoothing_penalty);
        __lhs_matrix *= std::pow(10.0, smoothing_penalty);

        __prepared = true;
    }

} //End FitPSF namespace.
