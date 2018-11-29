/**\file
 *
 * \brief Definitions of some methods of the PiecewiseBicubicMap class.
 */

#include "PiecewiseBicubicMap.h"

namespace PSF {

    const std::set<std::string> &
    PiecewiseBicubicMap::required_data_tree_quantities()
    {
        const std::string additional_required_data_tree_quantities[]={
#ifdef DEBUG
            "psffit.model",
#endif
            "psffit.grid",
            "psffit.psfmap"
        };

        static const std::set<std::string> required_data_tree_quantities=
            combine_required_tree_quantities(
                    additional_required_data_tree_quantities,
                    additional_required_data_tree_quantities
                    +
                    sizeof(additional_required_data_tree_quantities)
                    /
                    sizeof(additional_required_data_tree_quantities[0])
            );
        return required_data_tree_quantities;
    }

    template<class VECTOR_TYPE>
    void PiecewiseBicubicMap::fill_coefficients(
            const VECTOR_TYPE& coefficient_vector
    )
    {
        unsigned expansion_rows = ((__x_grid.size() - 2)
                                   *
                                   (__y_grid.size() - 2));

        if(expansion_rows == 0) {
            assert(coefficient_vector.size() == 0);
            set_num_terms(0);
            return;
        } else {
            set_num_terms(coefficient_vector.size() / 4 / expansion_rows);
        }
        unsigned matrix_step = coefficient_vector.size() / 4;
        assert(coefficient_vector.size() == 4 * matrix_step);

        const double *matrix_start=&(coefficient_vector[0]);
        __coefficients = matrix_start;
        __value_expansion = RowMajorMap(matrix_start,
                                        expansion_rows,
                                        num_terms());
        matrix_start += matrix_step;

        __x_deriv_expansion = RowMajorMap(matrix_start,
                                          expansion_rows,
                                          num_terms());
        matrix_start += matrix_step;

        __y_deriv_expansion = RowMajorMap(matrix_start,
                                          expansion_rows,
                                          num_terms());
        matrix_start += matrix_step;

        __xy_deriv_expansion = RowMajorMap(matrix_start,
                                           expansion_rows,
                                           num_terms());
        matrix_start += matrix_step;
    }

    PiecewiseBicubicMap::PiecewiseBicubicMap(const IO::H5IODataTree &data,
                                             double min_psf_span) :
        __coefficients(NULL)
    {
        std::cerr << "Reading model" << std::endl;
        if(
            data.get<std::string>("psffit.model",
                                  "",
                                  IO::translate_string) == "zero"
        ) {
            std::cerr << "model is 'zero'" << std::endl;
            __x_grid.resize(2);
            __y_grid.resize(2);
            __x_grid[0] = __y_grid[0] = -min_psf_span;
            __x_grid[1] = __y_grid[1] = min_psf_span;
        } else {
            std::cerr << "Reading grid" << std::endl;
            Grid grid = IO::parse_grid_string(
                data.get<std::string>("psffit.grid",
                                      "",
                                      IO::translate_string)
            );
            __x_grid=grid.x_grid;
            __y_grid=grid.y_grid;
        }

        std::cerr << "Reading coefficients" << std::endl;
        fill_coefficients(
            data.get<Eigen::VectorXd>(
                "psffit.psfmap",
                Eigen::VectorXd(),
                IO::TranslateToAny<Eigen::VectorXd>()
            )
        );
        if(__x_grid.front()>-min_psf_span) {
            __expanded_left=1;
            __x_grid.insert(__x_grid.begin(), -min_psf_span);
        } else __expanded_left=0;
        if(__x_grid.back()<min_psf_span) {
            __expanded_right=1;
            __x_grid.push_back(min_psf_span);
        } else __expanded_right=0;
        if(__y_grid.front()>-min_psf_span) {
            __expanded_down=1;
            __y_grid.insert(__y_grid.begin(), -min_psf_span);
        } else __expanded_down=0;
        if(__y_grid.back()<min_psf_span) {
            __expanded_up=1;
            __y_grid.push_back(min_psf_span);
        } else __expanded_up=0;
    }

    PiecewiseBicubicMap::PiecewiseBicubicMap(
            const Eigen::VectorXd &expansion_coef,
            const std::vector<double> &x_grid,
            const std::vector<double> &y_grid
    ) :
        __x_grid(x_grid),
        __y_grid(y_grid),
        __coefficients(NULL),
        __expanded_left(0),
        __expanded_right(0),
        __expanded_up(0),
        __expanded_down(0)
    {
        fill_coefficients(expansion_coef);
    }

    PiecewiseBicubic *PiecewiseBicubicMap::get_psf(
            const Eigen::VectorXd &terms,
            double background
    ) const
    {
        assert(terms.size() == num_terms());
        size_t
            unpadded_rows = (
                __y_grid.size() - 2 - __expanded_up - __expanded_down
            ),
            unpadded_cols = (
                __x_grid.size() - 2 - __expanded_left - __expanded_right
            ),
            padded_rows = __y_grid.size(),
            padded_cols = __x_grid.size(),
            padded_size = padded_rows * padded_cols;
        Eigen::VectorXd values = __value_expansion * terms,
                        x_deriv = __x_deriv_expansion * terms,
                        y_deriv = __y_deriv_expansion * terms,
                        xy_deriv = __xy_deriv_expansion * terms,
                        padded_values = Eigen::VectorXd::Constant(
                            padded_size,
                            background
                        ),
                        padded_x_deriv = Eigen::VectorXd::Zero(padded_size),
                        padded_y_deriv = Eigen::VectorXd::Zero(padded_size),
                        padded_xy_deriv = Eigen::VectorXd::Zero(padded_size);
        for(size_t row = 0; row < unpadded_rows; ++row) {
            size_t padded_start = (
                padded_cols * (row + 1 + __expanded_down)
                +
                1
                +
                __expanded_left
            ),
                   unpadded_start = unpadded_cols * row;
            padded_values.segment(padded_start, unpadded_cols).array() = (
                values.segment(unpadded_start, unpadded_cols).array()
                +
                background
            );
            padded_x_deriv.segment(padded_start, unpadded_cols) = (
                x_deriv.segment(unpadded_start, unpadded_cols)
            );
            padded_y_deriv.segment(padded_start, unpadded_cols) = (
                y_deriv.segment(unpadded_start, unpadded_cols)
            );
            padded_xy_deriv.segment(padded_start, unpadded_cols) = (
                xy_deriv.segment(unpadded_start, unpadded_cols)
            );
        }
        PiecewiseBicubic *result = new PiecewiseBicubic(
            __x_grid.begin(),
            __x_grid.end(),
            __y_grid.begin(),
            __y_grid.end()
        );
        result->set_values(padded_values.data(),
                           padded_x_deriv.data(),
                           padded_y_deriv.data(),
                           padded_xy_deriv.data());
        return result;
    }

} //End PSF namespace.
