/**\file
 *
 * \brief Define some of the methods of LinearSource class.
 *
 * \ingroup FitPSF
 */

#include "LinearSource.h"

namespace FitPSF {

    void LinearSource::calculate_predicted_pixel_values(
        double pixel_left,
        double pixel_bottom,
        const std::vector< std::valarray<double> > &parameter_sets,
        std::valarray<double> &pixel_integrals
    )
    {
        assert(16 * pixel_integrals.size() == parameter_sets[0].size());
        if (
            subpix_map().x_resolution() == 0
            &&
            subpix_map().y_resolution() == 0
        )
            pixel_integrals = __psf(pixel_left + 0.5,
                                    pixel_bottom + 0.5,
                                    parameter_sets);
        else {
#ifdef DEBUG
            assert(subpix_map().x_resolution() != 0);
            assert(subpix_map().y_resolution() != 0);
#endif
            double y_step = 1.0 / subpix_map().y_resolution(),
                   x_step = 1.0 / subpix_map().x_resolution();
            pixel_integrals = 0;
            double y0 = pixel_bottom + 0.5 / subpix_map().y_resolution();
            for(
                unsigned subpix_y = 0;
                subpix_y < subpix_map().y_resolution();
                ++subpix_y
            ) {
                double x0 = pixel_left + 0.5 / subpix_map().x_resolution();
                for(
                    unsigned subpix_x = 0;
                    subpix_x < subpix_map().x_resolution();
                    ++subpix_x
                ) {
                    if(subpix_map()(subpix_x, subpix_y))
                        pixel_integrals += (
                            subpix_map()(subpix_x, subpix_y)
                            *
                            __psf.integrate_rectangle_parameters(
                                x0, y0,
                                x_step, y_step,
                                parameter_sets
                            )
                        );
                    assert(!std::isnan(pixel_integrals.sum()));
                    x0 += x_step;
                }
                y0 += y_step;
            }
        }
    }

    void LinearSource::fill_pixel_integral_matrix(
        const std::vector< std::valarray<double> > &basis_parameter_sets,
        Eigen::Block<Eigen::MatrixXd> &shape_fit_integral_matrix,
        Eigen::VectorBlock<Eigen::VectorXd> &shape_fit_offset
    )
    {
        assert(shape_fit_integral_matrix.rows()
               ==
               static_cast<int>(shape_fit_pixel_count()));
        assert(shape_fit_offset.size()
               ==
               static_cast<int>(shape_fit_pixel_count()));

        assert(
            basis_parameter_sets[0].size()
            ==
            16 * static_cast<unsigned>(shape_fit_integral_matrix.cols() + 1)
        );

#ifdef VERBOSE_DEBUG
        std::cerr << "Basis parameter sets:" << std::endl;
        for(size_t set_i = 0; set_i < basis_parameter_sets.size(); ++set_i) {
            std::cerr << "\t";
            for(
                size_t par_i = 0;
                par_i < basis_parameter_sets[set_i].size();
                ++par_i
            )
                std::cerr << " " << basis_parameter_sets[set_i][par_i];
            std::cerr << std::endl;
        }
#endif

        unsigned num_param_sets = shape_fit_integral_matrix.cols() + 1;

        __flux_fit_integral_matrix.resize(
            flux_fit_pixel_count() - shape_fit_pixel_count(),
            num_param_sets - 1
        );
        __flux_fit_offset.resize(flux_fit_pixel_count()
                                 -
                                 shape_fit_pixel_count());

        unsigned shape_fit_index = 0,
                 flux_fit_index = 0;
        std::valarray<double> pixel_integrals(0.0, num_param_sets);
        for(
            ConstPixelIter pix_i = pixels().begin();
            pix_i != pixels().end();
            ++pix_i
        ) {
            calculate_predicted_pixel_values((*pix_i)->x() - x(),
                                             (*pix_i)->y() - y(),
                                             basis_parameter_sets,
                                             pixel_integrals);

#ifdef VERBOSE_DEBUG
            std::cerr << "Pixel(" 
                      << (*pix_i)->x()
                      << ", "
                      << (*pix_i)->y()
                      << ", shape: " 
                      << (*pix_i)->shape_fit()
                      << ", flux: "
                      << (*pix_i)->flux_fit()
                      << ") integrals:";
            for(size_t i = 0; i < pixel_integrals.size(); ++i)
                std::cerr << " " << pixel_integrals[i];
            std::cerr << "; var: " << (*pix_i)->variance() << std::endl;
#endif
            pixel_integrals /= std::sqrt((*pix_i)->variance());
            assert(!std::isnan(pixel_integrals.sum()));

            if((*pix_i)->shape_fit()) {
                assert((*pix_i)->flux_fit());

                shape_fit_integral_matrix.row(shape_fit_index) =
                    Eigen::Map<Eigen::RowVectorXd>(&(pixel_integrals[1]),
                                                   num_param_sets - 1);

                shape_fit_offset(shape_fit_index) = pixel_integrals[0];

                ++shape_fit_index;
            } else if((*pix_i)->flux_fit()) {
                __flux_fit_integral_matrix.row(flux_fit_index) =
                    Eigen::Map<Eigen::RowVectorXd>(&(pixel_integrals[1]),
                                                   num_param_sets - 1);

                __flux_fit_offset(flux_fit_index) = pixel_integrals[0];

                ++flux_fit_index;
            }
        }

        assert(shape_fit_index == shape_fit_pixel_count());
        assert(flux_fit_index
               ==
               flux_fit_pixel_count() - shape_fit_pixel_count());

        __shape_fit_integral_matrix =
            new Eigen::Block<Eigen::MatrixXd>(shape_fit_integral_matrix);
        __shape_fit_offset =
            new Eigen::VectorBlock<Eigen::VectorXd>(shape_fit_offset);
    }

    void LinearSource::fill_background_excess(
        Eigen::VectorBlock<Eigen::VectorXd> &shape_fit_background_excess
    )
    {
        unsigned excess_index = 0;
        for(
            ConstPixelIter pix_i = shape_fit_pixels_begin();
            pix_i != shape_fit_pixels_end();
            ++pix_i
        ) {
            shape_fit_background_excess[excess_index++] = background_excess(
                **pix_i,
                background_electrons(),
                background_electrons_variance()
            );
        }
        __shape_fit_bg_excess = new Eigen::VectorBlock<Eigen::VectorXd>(
            shape_fit_background_excess
        );

        __flux_fit_bg_excess.resize(flux_fit_pixel_count()
                                    -
                                    shape_fit_pixel_count());
        excess_index = 0;
        for(
            ConstPixelIter pix_i = flux_fit_pixels_begin();
            pix_i != flux_fit_pixels_end();
            ++pix_i
        ) {
            __flux_fit_bg_excess[excess_index++] = background_excess(
                **pix_i,
                background_electrons(),
                background_electrons_variance()
            );
        }
    }

    void LinearSource::pixel_excess_reductions(
        const Eigen::VectorXd &to_dot_with,
        double &dot_product,
        double &excess_sum_squares
    )
    {
        dot_product = (
            __shape_fit_bg_excess->dot(
                to_dot_with.head(shape_fit_pixel_count())
            )
            +
            __flux_fit_bg_excess.dot(
                to_dot_with.tail(__flux_fit_bg_excess.size())
            )
        );
        excess_sum_squares = (__shape_fit_bg_excess->squaredNorm()
                              +
                              __flux_fit_bg_excess.squaredNorm());
    }

    double LinearSource::fit_flux(const Eigen::VectorXd &psf_expansion_coef)
    {
        assert(overlaps().size() == 0);
        Eigen::VectorXd estimated_excesses(flux_fit_pixel_count());
        fill_fluxfit_column(
            psf_expansion_coef,
            estimated_excesses.head(shape_fit_pixel_count()),
            estimated_excesses.tail(flux_fit_pixel_count()
                                    -
                                    shape_fit_pixel_count())
        );
        return Source<PSF::PiecewiseBicubic>::fit_flux(
            estimated_excesses
        );
    }

    void LinearSource::prepare_for_fitting(
        const std::vector< std::valarray<double> > &basis_parameter_sets,
        Eigen::Block<Eigen::MatrixXd> shape_fit_integral_matrix,
        Eigen::VectorBlock<Eigen::VectorXd> shape_fit_offset,
        Eigen::VectorBlock<Eigen::VectorXd> shape_fit_background_excess
    ) 
    {
/*        OverlapSource<
            LinearSource,
            PSF::PiecewiseBicubic
        >::finalize_pixels();*/

        fill_pixel_integral_matrix(basis_parameter_sets,
                                   shape_fit_integral_matrix,
                                   shape_fit_offset);

        fill_background_excess(shape_fit_background_excess);

        __ready_to_fit = true;
    }

} //End FitPSF namespace.
