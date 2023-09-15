#include "PiecewiseBicubic.h"

namespace FitPSF {

    void time_this(const std::string
#ifdef TRACK_PROGRESS
                   message
#endif
                   = "")
    {
#ifdef TRACK_PROGRESS
        static clock_t previous_time;
        clock_t now = std::clock();
        if(message != "")
            std::cerr << message << " took "
                      << (now-previous_time) / CLOCKS_PER_SEC
                      << " sec." << std::endl;
        previous_time = now;
#endif
    }

    void fill_parameter_sets(
        const std::vector<double> &x_grid,
        const std::vector<double> &y_grid,
        std::vector< std::valarray<double> > &parameter_sets
    )
    {
        Core::vector_size_type x_res = x_grid.size() - 1,
                               y_res = y_grid.size() - 1,
                               num_sets = (x_res - 1) * (y_res - 1) * 4,
                               num_cells = x_res * y_res,
                               num_psf_param = (x_res + 1) * (y_res + 1);
        parameter_sets.resize(num_cells);
        for(Core::vector_size_type i = 0; i < num_cells; ++i)
            parameter_sets[i].resize(16 * num_sets, 0.0);
        std::valarray<double> psf_parameters(0.0, num_psf_param * 4);
        Core::vector_size_type psf_param_ind = 2 + x_res, set_ind = 0;
        for(
            Core::vector_size_type val_type_ind = 0;
            val_type_ind < 4;
            ++val_type_ind
        ) {
            for(Core::vector_size_type y_ind = 1; y_ind < y_res; ++y_ind) {
                double bot_height = y_grid[y_ind] - y_grid[y_ind - 1],
                       top_height = y_grid[y_ind + 1] - y_grid[y_ind];
                for(
                    Core::vector_size_type x_ind = 1;
                    x_ind < x_res;
                    ++x_ind
                ) {
                    psf_parameters[psf_param_ind] = 1;
                    Core::vector_size_type cell_psf_param_ind = (
                        x_ind - 1
                        +
                        (y_ind - 1) * (x_res + 1)
                    );
                    double left_width = x_grid[x_ind] - x_grid[x_ind - 1],
                           right_width = x_grid[x_ind + 1] - x_grid[x_ind];
                    PSF::calc_cell_coef(
                        &(psf_parameters[cell_psf_param_ind]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         2 * num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         3 * num_psf_param]),
                        x_res,
                        left_width,
                        bot_height,
                        parameter_sets[x_ind - 1 + (y_ind - 1) * x_res],
                        set_ind
                    );
                    cell_psf_param_ind += 1;
                    PSF::calc_cell_coef(
                        &(psf_parameters[cell_psf_param_ind]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         2 * num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         3 * num_psf_param]),
                        x_res,
                        right_width,
                        bot_height,
                        parameter_sets[x_ind + (y_ind - 1) * x_res],
                        set_ind
                    );
                    cell_psf_param_ind += x_res;
                    PSF::calc_cell_coef(
                        &(psf_parameters[cell_psf_param_ind]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         2 * num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         3 * num_psf_param]),
                        x_res,
                        left_width,
                        top_height,
                        parameter_sets[x_ind - 1 + y_ind * x_res],
                        set_ind
                    );
                    cell_psf_param_ind += 1;
                    PSF::calc_cell_coef(
                        &(psf_parameters[cell_psf_param_ind]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         2 * num_psf_param]),
                        &(psf_parameters[cell_psf_param_ind
                                         +
                                         3 * num_psf_param]),
                        x_res,
                        right_width,
                        top_height,
                        parameter_sets[x_ind + y_ind * x_res],
                        set_ind
                    );
                    psf_parameters[psf_param_ind] = 0;
                    set_ind += 16;
                    ++psf_param_ind;
                }
                psf_param_ind += 2;
            }
            psf_param_ind += 2 * (x_res + 1);
        }
    }

    void select_basis_vectors(
        const std::vector< std::valarray<double> > &parameter_sets,
        const PSF::PiecewiseBicubic &psf,
        Eigen::MatrixXd &basis
    )
    {
        unsigned num_param_sets = parameter_sets[0].size() / 16;
        Eigen::RowVectorXd full_integral_coef =
            Eigen::Map<Eigen::RowVectorXd>(
                &(
                    psf.integrate_rectangle_parameters(
                        (psf.min_x() + psf.max_x()) / 2.0,
                        (psf.min_y() + psf.max_y()) / 2.0,
                        psf.max_x() - psf.min_x(),
                        psf.max_y() - psf.min_y(),
                        parameter_sets
                    )[0]
                ),
                num_param_sets
            );
        Eigen::JacobiSVD< Eigen::MatrixXd,
                          Eigen::FullPivHouseholderQRPreconditioner >
            svd(full_integral_coef,
                Eigen::ComputeFullV | Eigen::ComputeFullU);
#ifdef DEBUG
        assert(svd.matrixU().size() == 1);
        assert(std::abs(svd.matrixU()(0)) == 1);
        assert(svd.singularValues().size() == 1);
        assert((svd.singularValues()(0) - full_integral_coef.norm())
               /
               full_integral_coef.norm()
               <
               10.0 * std::numeric_limits<double>::epsilon());
        assert(svd.matrixV().cols() == full_integral_coef.size());
        assert(
            (
                std::abs(full_integral_coef * svd.matrixV().col(0))
                -
                full_integral_coef.norm() * svd.matrixV().col(0).norm()
            )
            /
            std::abs(full_integral_coef * svd.matrixV().col(0))
            <
            10.0 * std::numeric_limits<double>::epsilon()
        );
        for(unsigned i = 1; i < svd.matrixV().cols(); ++i)
            assert(std::abs(full_integral_coef * svd.matrixV().col(i))
                   /
                   full_integral_coef.norm()
                   <
                   std::numeric_limits<double>::epsilon());
#endif
        if(svd.matrixU()(0, 0) < 0) basis = -svd.matrixV();
        else basis = svd.matrixV();
        basis.col(0) /= svd.singularValues()(0);
    }

    void fill_basis_parameter_sets(
        const std::vector< std::valarray<double> > &parameter_sets,
        Eigen::MatrixXd &basis,
        std::vector< std::valarray<double> > &basis_parameter_sets
    )
    {
#ifdef DEBUG
        assert(parameter_sets.size() == basis_parameter_sets.size());
        assert(basis.cols() * 16
               ==
               static_cast<int>(parameter_sets[0].size()));
        assert(basis.rows() * 16
               ==
               static_cast<int>(parameter_sets[0].size()));
#endif
        unsigned param_set_size = parameter_sets[0].size();
        for(unsigned cell = 0; cell < parameter_sets.size(); ++cell) {
            basis_parameter_sets[cell].resize(param_set_size, 0.0);
            for(unsigned col = 0; col < basis.cols(); ++col) {
                for(unsigned row = 0; row < basis.rows(); ++row) {
                    basis_parameter_sets[cell][std::slice(16 * col, 16, 1)]
                        +=
                        (
                            basis(row, col)
                            *
                            parameter_sets[cell][std::slice(16 * row, 16, 1)]
                        );
                }
            }
        }
    }

    void prepare_linear_regression(
        LinearSourceList& fit_sources,
        LinearSourceList& dropped_sources,
        const std::vector< std::valarray<double> > &basis_parameter_sets,
        Eigen::VectorXd &pixel_excesses,
        Eigen::MatrixXd &pixel_integral_matrix,
        Eigen::VectorXd &rhs_offset,
        Eigen::MatrixXd &symmetrized_pix_integral_matrix,
        Eigen::VectorXd &modified_pixel_excesses,
        Eigen::VectorXd &modified_rhs_offset
    )
    {
        unsigned total_pix_ind = 0,
                 source_ind = 0,
                 num_param_sets = basis_parameter_sets[0].size() / 16;

        assert(symmetrized_pix_integral_matrix.cols() == num_param_sets - 1);
        assert(symmetrized_pix_integral_matrix.rows() ==
                static_cast<int>((num_param_sets - 1) * fit_sources.size()));
        assert(modified_pixel_excesses.size() ==
                static_cast<int>((num_param_sets - 1) * fit_sources.size()));
        assert(pixel_integral_matrix.cols() ==
                static_cast<int>(num_param_sets - 1));
        assert(modified_rhs_offset.size() == modified_pixel_excesses.size());
        assert(rhs_offset.size() == pixel_excesses.size());


        std::valarray<double> pixel_integrals(0.0, num_param_sets);
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            unsigned source_pixels = (*si)->shape_fit_pixel_count();
#ifdef VERBOSE_DEBUG
            std::cerr << "Preparing source at ("
                      << (*si)->x() << ", " << (*si)->y()
                      << "), containing " << source_pixels
                      << " shape fitting out of  "
                      << (*si)->flux_fit_pixel_count()
                      << " pixels for fitting!"
                      << std::endl;
#endif
            (*si)->prepare_for_fitting(
                basis_parameter_sets,
                pixel_integral_matrix.block(total_pix_ind,
                                            0,
                                            source_pixels,
                                            num_param_sets -1),
                rhs_offset.segment(total_pix_ind, source_pixels),
                pixel_excesses.segment(total_pix_ind, source_pixels)
            );

            symmetrized_pix_integral_matrix.block(
                source_ind * (num_param_sets - 1),
                0,
                num_param_sets - 1,
                num_param_sets - 1
            ) = (
                pixel_integral_matrix.block(total_pix_ind,
                                            0,
                                            source_pixels,
                                            num_param_sets - 1).transpose()
                *
                pixel_integral_matrix.block(total_pix_ind,
                                            0,
                                            source_pixels,
                                            num_param_sets - 1)
            );

            modified_pixel_excesses.segment(
                source_ind * (num_param_sets - 1),
                num_param_sets - 1
            ) = (
                pixel_integral_matrix.block(total_pix_ind,
                                            0,
                                            source_pixels,
                                            num_param_sets - 1).transpose()
                *
                pixel_excesses.segment(total_pix_ind, source_pixels)
            );

            modified_rhs_offset.segment(source_ind * (num_param_sets - 1),
                                        num_param_sets - 1)
                =
                pixel_integral_matrix.block(total_pix_ind,
                                            0,
                                            source_pixels,
                                            num_param_sets - 1).transpose()
                *
                rhs_offset.segment(total_pix_ind, source_pixels);
            total_pix_ind += source_pixels;

            ++source_ind;
        }
        assert(total_pix_ind == pixel_integral_matrix.rows());

        for(
            LinearSourceList::iterator si = dropped_sources.begin();
            si != dropped_sources.end();
            ++si
        ) {
            (*si)->finalize_pixels();

            (*si)->prepare_for_fitting(
                basis_parameter_sets,
                pixel_integral_matrix.bottomRows(0),
                rhs_offset.tail(0),
                pixel_excesses.tail(0)
            );
        }

        assert(!std::isnan(pixel_integral_matrix.sum()));
        assert(!std::isnan(rhs_offset.sum()));
        assert(!std::isnan(symmetrized_pix_integral_matrix.sum()));
        assert(!std::isnan(modified_pixel_excesses.sum()));
        assert(!std::isnan(modified_rhs_offset.sum()));
    }

    void fill_poly_coef_matrix(const LinearSourceList &fit_sources,
                               Eigen::MatrixXd &poly_coef_matrix)
    {
        unsigned
            source_ind = 0,
            num_poly_terms = fit_sources.front()->expansion_terms().size();
#ifdef DEBUG
        assert(num_poly_terms == poly_coef_matrix.cols());
        assert(static_cast<int>(num_poly_terms * fit_sources.size())
               ==
               poly_coef_matrix.rows());
#endif

        for(
            LinearSourceList::const_iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            poly_coef_matrix.block(
                source_ind * num_poly_terms,
                0,
                num_poly_terms,
                num_poly_terms
            ) = (
                (*si)->expansion_terms()
                *
                (*si)->expansion_terms().transpose()
            );
            ++source_ind;
        }
    }

    void fill_matrix_to_invert(
        const LinearSourceList &fit_sources,
        const Eigen::MatrixXd &symmetrized_pix_integral_matrix,
        const Eigen::MatrixXd &poly_coef_matrix,
        Eigen::MatrixXd &matrix_to_invert
    )
    {
        unsigned num_psf_params = symmetrized_pix_integral_matrix.cols(),
                 num_poly_terms = poly_coef_matrix.cols();
#ifdef DEBUG
        assert(symmetrized_pix_integral_matrix.rows()
               ==
               static_cast<int>(fit_sources.size() * num_psf_params));
        assert(matrix_to_invert.rows()
               ==
               static_cast<int>(num_psf_params * num_poly_terms));
        assert(matrix_to_invert.rows() == matrix_to_invert.cols());
#endif

        for(unsigned row = 0; row < num_psf_params; ++row)
            for(unsigned col = 0; col < num_psf_params; ++col) {
                unsigned source_ind = 0;
                for(
                    LinearSourceList::const_iterator
                        si = fit_sources.begin();
                    si != fit_sources.end();
                    ++si
                ) {
                    matrix_to_invert.block(
                        row * num_poly_terms,
                        col * num_poly_terms,
                        num_poly_terms,
                        num_poly_terms
                    ) += (
                        std::pow((*si)->flux(0).value(), 2)
                        *
                        symmetrized_pix_integral_matrix(
                            source_ind * num_psf_params + row,
                            col
                        )
                        *
                        poly_coef_matrix.block(
                            source_ind * num_poly_terms,
                            0,
                            num_poly_terms,
                            num_poly_terms
                        )
                    );
                    ++source_ind;
                }
            }
    }

    void fill_flux_scaled_modified_rhs(
        const LinearSourceList &fit_sources,
        const Eigen::VectorXd &modified_pixel_excesses,
        const Eigen::VectorXd &modified_rhs_offset,
        Eigen::VectorXd &flux_scaled_modified_rhs
    )
    {
        unsigned
            num_poly_terms = fit_sources.front()->expansion_terms().size(),
            num_param_sets = (modified_pixel_excesses.size()
                              /
                              fit_sources.size()),
            source_ind = 0;
#ifdef DEBUG
        assert(modified_pixel_excesses.size() % fit_sources.size() == 0);
        assert(flux_scaled_modified_rhs.size()
               ==
               num_param_sets * num_poly_terms);
#endif

        for(
            LinearSourceList::const_iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
#ifdef DEBUG
            assert((*si)->flux(0).value() > 0);
#endif
            for(
                unsigned param_set_ind = 0;
                param_set_ind < num_param_sets;
                ++param_set_ind
            ) {
                flux_scaled_modified_rhs.segment(
                    num_poly_terms * param_set_ind,
                    num_poly_terms
                ) += (
                    (*si)->flux(0).value()
                    *
                    (
                        modified_pixel_excesses(num_param_sets * source_ind
                                                +
                                                param_set_ind)
                        -
                        (*si)->flux(0).value()
                        *modified_rhs_offset(num_param_sets * source_ind
                                             +
                                             param_set_ind)
                    )
                    *
                    (*si)->expansion_terms()
                );
            }
            ++source_ind;
        }
    }

    void estimate_initial_amplitudes(LinearSourceList &fit_sources,
                                     double gain)
    {
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
#ifdef VERBOSE_DEBUG
            std::cerr << "Source ("
                      << *si
                      << ") at ("
                      << (*si)->x()
                      << ", "
                      << (*si)->y()
                      << ") background = "
                      << (*si)->background_electrons()
                      << ", # flux fit pixels: "
                      << (*si)->flux_fit_pixel_count()
                      << std::endl;
#endif
            double total_flux_electrons = 0, total_var_electrons2 = 0;
            for(
                ConstPixelIter pix_i = (*si)->pixels().begin();
                pix_i != (*si)->pixels().end();
                ++pix_i
            ) {
                total_flux_electrons += ((*pix_i)->measured()
                                         -
                                         (*si)->background_electrons());
                total_var_electrons2 += (*pix_i)->variance();
#ifdef VERBOSE_DEBUG
                std::cerr << "With pix ("
                          << *pix_i
                          << ") at ("
                          << (*pix_i)->x()
                          << ", "
                          << (*pix_i)->y()
                          << "), "
                          << (*pix_i)->measured()
                          << " flux = "
                          << total_flux_electrons
                          << ", flux fit: "
                          << (*pix_i)->flux_fit()
                          << std::endl;
#endif
            }
            (*si)->flux(0).value() = total_flux_electrons / gain;
            (*si)->flux(0).error() = std::sqrt(total_var_electrons2) / gain;
            (*si)->flux(0).flag() = (*si)->quality_flag();
            assert(
                (*si)->flux_fit_pixel_count() == 0
                ||
                !std::isnan((*si)->flux(0).value())
            );
        }
    }

    ///Use the given PSF map to estimate initial fluxes for the sources.
    void estimate_initial_amplitudes(LinearSourceList &fit_sources,
                                     OverlapGroupList &overlap_groups,
                                     const PSF::PiecewiseBicubicMap &psf_map)
    {
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            if((*si)->overlaps().size()) continue;
            PSF::PiecewiseBicubic
                *psf = psf_map.get_psf((*si)->expansion_terms());
            (*si)->fit_flux(*psf);
            delete psf;
        }

        for(
            OverlapGroupList::iterator group_i = overlap_groups.begin();
            group_i != overlap_groups.end();
            ++group_i
        )
            group_i->fit_fluxes(psf_map);
    }

    void estimate_initial_amplitudes(
        LinearSourceList &fit_sources,
        const Core::SubPixelMap &subpix_map,
        const Core::Image<double> &observed_image,
        double gain,
        double aperture
    )
    {
        double more_aperture = 1.5 * aperture;
        std::vector<double> grid(2),
                            zeros(4, 0),
                            psf_values(4, 0.25 / std::pow(more_aperture, 2));
        grid[0] = -more_aperture;
        grid[1] = more_aperture;
        PSF::PiecewiseBicubic psf(grid.begin(),
                                  grid.end(),
                                  grid.begin(),
                                  grid.end());
        psf.set_values(psf_values.begin(),
                       zeros.begin(),
                       zeros.begin(),
                       zeros.begin());
        Core::SubPixelMap uniform_map(1, 1, "");
        uniform_map(0, 0) = 1;
        Core::SubPixelCorrectedFlux<Core::SubPixelMap> measure_flux(
            observed_image,
            (subpix_map.x_resolution() == 0 ? uniform_map : subpix_map),
            0,
            std::list<double>(1, aperture),
            gain
        );
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            (*si)->flux(0) = measure_flux(
                (*si)->x(),
                (*si)->y(),
                psf,
                (*si)->background_electrons() / gain,
                std::sqrt((*si)->background_electrons_variance()) / gain
            )[0];
            (*si)->flux(0).flag() = (*si)->quality_flag();
#ifdef DEBUG
            assert(!std::isnan((*si)->flux(0).value()));
#endif
        }
    }

    /*
    void estimate_initial_amplitudes(LinearSourceList &fit_sources,
            const Core::SubPixelMap &subpix_map,
            const Core::Image<double> &observed_image,
            double gain, double aperture,
            unsigned poly_order, double x_offset, double y_offset,
            double x_scale, double y_scale, double max_abs_amplitude_change,
            double max_rel_amplitude_change)
    {
        std::vector<double> grid(3);
        grid[0]=-aperture;
        grid[1]=0;
        grid[2]=aperture;
        PSF::PiecewiseBicubic dummy_psf(grid.begin(), grid.end(),
                                        grid.begin(), grid.end());
        estimate_initial_amplitudes(fit_sources, subpix_map, observed_image,
                gain, aperture);
        std::vector< std::valarray<double> > parameter_sets;
        fill_parameter_sets(grid, grid, parameter_sets);
        std::vector<double> zeros(9, 0);
        dummy_psf.set_values(zeros.begin(), zeros.begin(), zeros.begin(),
                zeros.begin());
        size_t max_source_pixels,
               num_pixels=count_pixels(fit_sources, &max_source_pixels);
        size_t num_psf_params=4,
               num_poly_terms=num_poly_coef(poly_order);
        Eigen::VectorXd rhs(num_pixels),
            modified_pixel_excesses(fit_sources.size()*num_psf_params),
            modified_rhs_offset;
        fill_pixel_excesses(fit_sources, pixel_excesses);
        Eigen::MatrixXd
            pixel_integral_matrix(num_pixels, num_psf_params),
            symmetrized_pix_integral_matrix(num_psf_params*fit_sources.size(),
                    num_psf_params),
            poly_coef_matrix(fit_sources.size()*num_poly_terms,
                    num_poly_terms);
        prepare_linear_regression(fit_sources, subpix_map, parameter_sets,
                dummy_psf, rhs, pixel_integral_matrix,
                symmetrized_pix_integral_matrix, modified_pixel_excesses,
                modified_rhs_offset);
        time_this("Preparing initial amplitude the linear regression");
        fill_poly_coef_matrix(fit_sources, poly_order, x_offset, y_offset,
                x_scale, y_scale, poly_coef_matrix);
        time_this("Filling the initial amplitude polynomial coefficient matrix");
        double amplitude_change=1, max_amplitude_change=0;
        Eigen::VectorXd best_fit_poly_coef;
        while(amplitude_change>max_amplitude_change)
        {
            amplitude_change=fit_piecewise_bicubic_psf_step(fit_sources,
                    symmetrized_pix_integral_matrix, poly_coef_matrix,
                    full_psf_integral_matrix, pixel_integral_matrix, rhs,
                    modified_pixel_excesses, poly_order, x_offset, y_offset,
                    x_scale, y_scale, best_fit_poly_coef);
            double amplitude_vector_norm=0;
            for(LinearSourceList::iterator
                    si=fit_sources.begin(); si!=fit_sources.end(); ++si)
                amplitude_vector_norm+=std::pow((*si)->flux(0).value(),2);
            amplitude_vector_norm=std::sqrt(amplitude_vector_norm);
            max_amplitude_change=std::max(max_abs_amplitude_change,
                    max_rel_amplitude_change*amplitude_vector_norm);
#ifdef VERBOSE_DEBUG
            std::cerr << "Amplitude vector norm=" << amplitude_vector_norm
                << std::endl;
#endif
        };
    }*/

    double update_fluxes(LinearSourceList &fit_sources,
                         OverlapGroupList &overlap_groups,
                         const Eigen::VectorXd &best_fit)
    {
        double result = 0;
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        )
            if((*si)->overlaps().size() == 0) {
                result += std::pow((*si)->fit_flux(best_fit), 2);
#ifdef VERBOSE_DEBUG
                if(std::isnan(result))
                    std::cout << "NaN flux for Source(x = " << (*si)->x()
                              << ", y = " << (*si)->y()
                              << "): flux=" << (*si)->flux(0).value()
                              << ", shape fit pixel count="
                              << (*si)->shape_fit_pixel_count()
                              << ", flux fit pixel count="
                              << (*si)->flux_fit_pixel_count()
                              << ", chi2="
                              << (*si)->chi2()
                              << std::endl;
#endif
                assert(!std::isnan(result));
            }
        for(
            OverlapGroupList::iterator gi = overlap_groups.begin();
            gi != overlap_groups.end();
            ++gi
        ) {
            result += gi->fit_fluxes(best_fit);
        }

        return std::sqrt(result);
    }

    double fit_piecewise_bicubic_psf_step(
        LinearSourceList &fit_sources,
        OverlapGroupList &overlap_groups,
        const Eigen::MatrixXd &symmetrized_pix_integral_matrix,
        const Eigen::MatrixXd &poly_coef_matrix,
        const Eigen::VectorXd &modified_pixel_excesses,
        const Eigen::VectorXd &modified_rhs_offset,
        const PiecewiseBicubicPSFSmoothing &smoothing,
        Eigen::VectorXd &best_fit
    )
    {
        unsigned num_psf_params = symmetrized_pix_integral_matrix.cols(),
                 num_fit_params = poly_coef_matrix.cols() * num_psf_params;
        Eigen::MatrixXd matrix_to_invert =
            Eigen::MatrixXd::Zero(num_fit_params, num_fit_params);
        Eigen::VectorXd flux_scaled_modified_rhs =
            Eigen::VectorXd::Zero(num_fit_params);

        time_this();
        fill_matrix_to_invert(fit_sources,
                              symmetrized_pix_integral_matrix,
                              poly_coef_matrix,
                              matrix_to_invert);
        assert(!std::isnan(matrix_to_invert.sum()));
        time_this("Creating matrix to invert");

        fill_flux_scaled_modified_rhs(fit_sources,
                                      modified_pixel_excesses,
                                      modified_rhs_offset,
                                      flux_scaled_modified_rhs);
        assert(!std::isnan(flux_scaled_modified_rhs.sum()));
        time_this("Filling the final RHS");

        if(smoothing) {
#ifdef VERBOSE_DEBUG
            std::cerr << "Smoothing matrix size: "
                      << smoothing.lhs_correction().rows()
                      << "x"
                      << smoothing.lhs_correction().cols()
                      << std::endl;
            std::cerr << "fitting matrix size: "
                      << matrix_to_invert.rows()
                      << "x"
                      << matrix_to_invert.cols()
                      << std::endl;
#endif
            matrix_to_invert += smoothing.lhs_correction();
            flux_scaled_modified_rhs += smoothing.rhs_correction();
            time_this("Adding smoothing constraint");
        }

        Eigen::LDLT<Eigen::MatrixXd> decomposition(matrix_to_invert);
        time_this("Deriving the decomposition");

        best_fit = decomposition.solve(flux_scaled_modified_rhs);
        assert(!std::isnan(best_fit.sum()));
#ifdef VERBOSE_DEBUG
        std::cerr << "Best fit: " << best_fit << std::endl;
#endif

        time_this("Applying the decomposition");

        double result = update_fluxes(fit_sources, overlap_groups, best_fit);
        time_this("Updating the fluxes");
#ifdef TRACK_PROGRESS
        std::cerr << "Amplitude change=" << result << std::endl;
#endif
        return result;
    }

    size_t count_pixels(const LinearSourceList &fit_sources,
                        size_t *max_source_pixels)
    {
        if(max_source_pixels) *max_source_pixels = 0;
        size_t num_pixels = 0;
        for(
            LinearSourceList::const_iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            (*si)->finalize_pixels();
            num_pixels += (*si)->shape_fit_pixel_count();
            if(max_source_pixels)
                *max_source_pixels=std::max((*si)->shape_fit_pixel_count(),
                                            *max_source_pixels);
        }
        return num_pixels;
    }

    ///\brief Discards pixels with too large residuals after the last PSF
    ///fit.
    ///
    ///Returns the number of pixels discarded.
    unsigned discard_outlier_pixels(
        ///The sources participating in the fit.
        LinearSourceList &fit_sources,

        ///See prepare_linear_regression
        const Eigen::MatrixXd &pixel_integral_matrix,

        ///See fill_fit_pixel_excesses.
        const Eigen::VectorXd &pixel_excesses,

        ///See prepare_linear_regression.
        const Eigen::VectorXd &rhs_offset,

        ///If source pixels with residuals larger than this are found, they
        ///are discarded from the fit.
        double max_residual,

        ///The current best fit polynomial expansion coefficients for the PSF
        ///parameters.
        const Eigen::VectorXd &best_fit)
    {
        unsigned num_psf_params = pixel_integral_matrix.cols(),
                 pix_ind = 0,
                 result = 0;
#ifdef DEBUG
        assert(pixel_excesses.size() == rhs_offset.size());
        assert(pixel_integral_matrix.rows()
               ==
               static_cast<int>(pixel_excesses.size()));
#endif

        Eigen::VectorXd psf_params(num_psf_params);
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            unsigned shape_fit_pixel_count = (*si)->shape_fit_pixel_count();
            (*si)->fill_psf_params(best_fit, psf_params);
            Eigen::VectorXd predicted_excesses =
                (*si)->flux(0).value()
                *
                (
                    pixel_integral_matrix.block(pix_ind,
                                                0,
                                                shape_fit_pixel_count,
                                                num_psf_params)
                    *
                    psf_params
                    +
                    rhs_offset.segment(pix_ind, shape_fit_pixel_count)
                ),
                abs_residuals = (pixel_excesses.segment(pix_ind,
                                                        shape_fit_pixel_count)
                                 -
                                 predicted_excesses).cwiseAbs();
            if(abs_residuals.maxCoeff() > max_residual) {
                int residual_ind = 0;
                for(
                    PixelIter pix_i = (*si)->pixels().begin();
                    pix_i != (*si)->pixels().end();
                    ++pix_i
                ) {
                    if((*pix_i)->shape_fit()) {
                        assert((*pix_i)->flux_fit());
                        if(abs_residuals[residual_ind] > max_residual) {
                            (*pix_i)->exclude_from_shape_fit();
                            (*pix_i)->exclude_from_flux_fit();
                            ++result;
                        }
                        ++residual_ind;
                    }
                }
                assert(residual_ind == abs_residuals.size());
            }
            pix_ind += shape_fit_pixel_count;
        }
        return result;
    }

    ///\brief Discards sources unsuitable for fitting.
    ///
    ///At least one of the following is true of the discarded sources:
    /// - the flux estimate is nan or negative
    /// - the chi squared is too big (only if discard_by_chi2 is true).
    /// - less than two pixels are assigned to the source
    unsigned discard_sources(
        ///The sources participating in the fit.
        LinearSourceList &fit_sources,

        ///The maximum allowed reduced chi squared before discarding.
        double max_chi2,

        ///Should sources be discarded besed on their chi squared.
        bool discard_by_chi2,

        ///The sources dropped from the fit.
        LinearSourceList &dropped_sources
    )
    {
        size_t num_discarded_sources = 0;
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
        ) {
            if(
                std::isnan((*si)->flux(0).value())
                ||
                (*si)->flux(0).value() < 0
                ||
                (*si)->shape_fit_pixel_count() < 2
                ||
                (
                    discard_by_chi2
                    &&
                    (
                        (*si)->chi2()
                        /
                        ((*si)->flux_fit_pixel_count() - 1)
                    ) > max_chi2
                )
            ) {
#ifdef TRACK_PROGRESS
                std::cout << "Discarding source(x = " << (*si)->x()
                          << ", y = " << (*si)->y()
                          << "): flux=" << (*si)->flux(0).value()
                          << ", shape fit pixel count="
                          << (*si)->shape_fit_pixel_count()
                          << ", flux fit pixel count="
                          << (*si)->flux_fit_pixel_count()
                          << ", chi2="
                          << (*si)->chi2()
                          << std::endl;
#endif
                LinearSourceList::iterator to_drop = si++;
                (*to_drop)->exclude_from_shape_fit();
                dropped_sources.splice(dropped_sources.end(),
                                       fit_sources,
                                       to_drop);
                ++num_discarded_sources;
            } else ++si;
        }
        return num_discarded_sources;
    }

    ///\brief Discards pixels whose overlap with the PSF domain has zero
    ///sub-pixel sensitivity.
    ///
    ///Returns the number of pixels discarded.
    unsigned discard_zero_flux_pixels(
        ///The sources participating in the fit.
        LinearSourceList &fit_sources,

        ///The sub-pixel sensitivity map.
        const Core::SubPixelMap &subpix_map,

        ///A PSF with the correct grid of cells.
        const PSF::PiecewiseBicubic &psf
    )
    {
        if(subpix_map.x_resolution() != 0 && subpix_map.min() > 0) return 0;
#ifdef DEBUG
        assert(
                (
                    subpix_map.x_resolution() == 0
                    &&
                    subpix_map.y_resolution() == 0
                )
                ||
                subpix_map.min() == 0
        );
#endif
        unsigned result = 0;
        double y_step = 1.0 / subpix_map.y_resolution(),
               x_step = 1.0 / subpix_map.x_resolution();
        std::vector<double> x_grid(2), y_grid(2), ones(4, 1.0), zeros(4, 0);
        x_grid[0] = psf.min_x();
        x_grid[1] = psf.max_x();
        y_grid[0] = psf.min_y();
        y_grid[1] = psf.max_y();
        PSF::PiecewiseBicubic unit_psf(x_grid.begin(),
                                       x_grid.end(),
                                       y_grid.begin(),
                                       y_grid.end());
        unit_psf.set_values(ones.begin(), zeros.begin(), zeros.begin(),
                            zeros.begin());
        for(
            LinearSourceList::iterator si = fit_sources.begin();
            si != fit_sources.end();
            ++si
        ) {
            for(
                PixelIter pix_i = (*si)->pixels().begin();
                pix_i != (*si)->pixels().end();
                ++pix_i
            ) {
                if(!(*pix_i)->shape_fit()) continue;
                assert((*pix_i)->flux_fit());
                bool nonzero = false;
                if(
                    subpix_map.x_resolution() == 0
                    &&
                    subpix_map.y_resolution() == 0
                ) {
                    double x0 = (*pix_i)->x() + 0.5 - (*si)->x(),
                           y0 = (*pix_i)->y() + 0.5 - (*si)->y();
                    nonzero = (x0 > psf.min_x() && x0 < psf.max_x()
                               &&
                               y0 > psf.min_y() && y0 < psf.max_y());
                } else {
                    double y0 = ((*pix_i)->y()
                                 -
                                 (*si)->y()
                                 +
                                 0.5 / subpix_map.y_resolution());
                    for(
                        unsigned subpix_y = 0;
                        subpix_y < subpix_map.y_resolution();
                        ++subpix_y
                    ) {
                        double x0 = ((*pix_i)->x()
                                     -
                                     (*si)->x()
                                     +
                                     0.5 / subpix_map.x_resolution());
                        for(
                            unsigned subpix_x = 0;
                            subpix_x < subpix_map.x_resolution();
                            ++subpix_x
                        ) {
#ifdef VERBOSE_DEBUG
                            std::cerr << "subpix("
                                      << subpix_x
                                      << ", "
                                      << subpix_y
                                      << ") = "
                                      << subpix_map(subpix_x, subpix_y)
                                      << "; "
                                      << "PSF integral ("
                                      << x0 - x_step / 2.0
                                      << " < x < "
                                      << x0 + x_step / 2.0
                                      << ", "
                                      << y0 - y_step / 2.0
                                      << " < x < "
                                      << y0 + y_step / 2.0
                                      << ") = "
                                      << unit_psf.integrate(x0,
                                                            y0,
                                                            x_step,
                                                            y_step)
                                      << std::endl;
#endif
                            if(
                                subpix_map(subpix_x, subpix_y) != 0
                                &&
                                unit_psf.integrate(x0,
                                                   y0,
                                                   x_step,
                                                   y_step) != 0
                            ) {
                                nonzero = true;
                                break;
                            }
                            x0 += x_step;
                        }
                        if(nonzero) break;
                        y0 += y_step;
                    }
                }
#ifdef VERBOSE_DEBUG
                std::cerr << "Pixel("
                          << (*pix_i)->x()
                          << ", "
                          << (*pix_i)->y()
                          << ") is nonzero: "
                          << nonzero
                          << std::endl;
#endif
                if(!nonzero) {
                    (*pix_i)->exclude_from_shape_fit();
                    (*pix_i)->exclude_from_flux_fit();
                    ++result;
                }
            }
        }
        return result;
    }

    ///Transforms a set of best fit coefficients from some basis to directly
    ///expressing the spatial dependence of PSF parameters.
    void transform_best_fit_to_direct_parameters(
        ///The matrix of basis vectors used when solving
        const Eigen::MatrixXd &basis,

        ///The best fit coefficients in the above basis. Overwritten on
        ///output!
        Eigen::VectorXd &best_fit_poly_coef
    )
    {
        typedef Eigen::Map< Eigen::VectorXd,
                            Eigen::Unaligned,
                            Eigen::InnerStride<>
                          > SingleTermMap;
        unsigned num_spatial_terms = best_fit_poly_coef.size()
                                     /
                                     (basis.cols() - 1);
        Eigen::VectorXd result(best_fit_poly_coef.size() + num_spatial_terms);
        for(unsigned offset = 0; offset < num_spatial_terms; ++offset) {
            SingleTermMap basis_single_spatial_term_coef(
                                  &best_fit_poly_coef(offset),
                                  basis.cols() - 1,
                                  Eigen::InnerStride<>(num_spatial_terms)),
                          result_single_spatial_term_coef(
                                  &result(offset),
                                  basis.cols(),
                                  Eigen::InnerStride<>(num_spatial_terms));
            result_single_spatial_term_coef = (
                basis.rightCols(basis.cols() - 1)
                *
                basis_single_spatial_term_coef
            );
            if(offset == 0) result_single_spatial_term_coef += basis.col(0);
        }
        best_fit_poly_coef = result;
    }

    void convert_flux_to_adu( LinearSourceList& sources, double gain )
    {
        for (
            LinearSourceList::iterator s_i = sources.begin();
            s_i != sources.end();
            ++s_i
        ) {
            Core::Flux& flux = (*s_i)->flux(0);

            flux.value() /= gain;
            flux.error() /= gain;
        }
    }

    ///\brief PSF fitting, but assumed that sources have their initial
    ///amplitudes set.
    bool fit_piecewise_bicubic_psf(LinearSourceList &fit_sources,
                                   OverlapGroupList &overlap_groups,
                                   LinearSourceList &dropped_sources,
                                   const std::vector<double> &x_grid,
                                   const std::vector<double> &y_grid,
                                   const Core::SubPixelMap &subpix_map,
                                   double max_abs_amplitude_change,
                                   double max_rel_amplitude_change,
                                   double max_chi2,
                                   double pixel_rejection,
                                   double min_convergence_rate,
                                   int max_iterations,
                                   double smoothing_penalty,
                                   Eigen::VectorXd &best_fit_poly_coef,
                                   double gain)
    {
        time_this();
        PSF::PiecewiseBicubic dummy_psf(x_grid.begin(),
                                        x_grid.end(),
                                        y_grid.begin(),
                                        y_grid.end());
        std::vector<double> zeros(x_grid.size() * y_grid.size(), 0);
        dummy_psf.set_values(zeros.begin(), zeros.begin(),
                             zeros.begin(), zeros.begin());

        size_t num_discarded_sources = 0,
               discarded_pixels = discard_zero_flux_pixels(fit_sources,
                                                           subpix_map,
                                                           dummy_psf);

        std::vector< std::valarray<double> > parameter_sets;
        fill_parameter_sets(x_grid, y_grid, parameter_sets);
        time_this("Filling parameter sets");
        Eigen::MatrixXd basis;
        select_basis_vectors(parameter_sets, dummy_psf, basis);
        time_this("Selecting basis vectors.");

        std::vector< std::valarray<double> >
            basis_parameter_sets(parameter_sets.size());
        fill_basis_parameter_sets(parameter_sets,
                                  basis,
                                  basis_parameter_sets);
        time_this("Transforming parameter sets.");

        PiecewiseBicubicPSFSmoothing smoothing;
        if(std::isfinite(smoothing_penalty)) {
            smoothing.prepare_smoothing(smoothing_penalty,
                                        fit_sources,
                                        x_grid,
                                        y_grid,
                                        basis_parameter_sets);
            time_this("Preparing smoothing");
        } else
            time_this("No smoothing");

        bool discard_by_chi2 = false;
        bool converged = true;
        discard_sources(fit_sources,
                        max_chi2,
                        discard_by_chi2,
                        dropped_sources);
        time_this("Initial discarding of sources");
        while(
            num_discarded_sources > 0
            || discarded_pixels > 0
            || !discard_by_chi2
        ) {
            size_t max_source_pixels,
                   num_pixels = count_pixels(fit_sources,
                                             &max_source_pixels);
            double max_residual = (pixel_rejection>0
                                   ? pixel_rejection
                                   : Core::Inf);
#ifdef TRACK_PROGRESS
            std::cerr << "After discarding: " << discarded_pixels
                      << " pixels and " << num_discarded_sources
                      << " sources, " << fit_sources.size()
                      << " sources remain:" << std::endl;
#ifndef NDEBUG
            for(
                LinearSourceList::iterator src_i = fit_sources.begin();
                src_i != fit_sources.end();
                ++src_i
            )
                std::cerr << "\t" << (*src_i)->id() << std::endl;
#endif
#endif
            size_t num_psf_params = ((x_grid.size() - 2)
                                     *
                                     (y_grid.size() - 2) * 4
                                     -
                                     1),
                   num_poly_terms = (
                       fit_sources.front()->expansion_terms().size()
                   );
            Eigen::VectorXd
                pixel_excesses(num_pixels),
                rhs_offset(num_pixels),
                modified_pixel_excesses(fit_sources.size() * num_psf_params),
                modified_rhs_offset(modified_pixel_excesses.size());
            Eigen::MatrixXd
                pixel_integral_matrix(num_pixels, num_psf_params),
                symmetrized_pix_integral_matrix(
                        num_psf_params * fit_sources.size(),
                        num_psf_params
                ),
                poly_coef_matrix(
                        fit_sources.size() * num_poly_terms,
                        num_poly_terms
                );
#ifndef NDEBUG
            pixel_excesses.setConstant(Core::NaN);
            rhs_offset.setConstant(Core::NaN);
            modified_pixel_excesses.setConstant(Core::NaN);
            modified_rhs_offset.setConstant(Core::NaN);
            pixel_integral_matrix.setConstant(Core::NaN);
            symmetrized_pix_integral_matrix.setConstant(Core::NaN);
            poly_coef_matrix.setConstant(Core::NaN);
#endif
            prepare_linear_regression(fit_sources,
                                      dropped_sources,
                                      basis_parameter_sets,
                                      pixel_excesses,
                                      pixel_integral_matrix,
                                      rhs_offset,
                                      symmetrized_pix_integral_matrix,
                                      modified_pixel_excesses,
                                      modified_rhs_offset);
            time_this("Preparing the linear regression");
            fill_poly_coef_matrix(fit_sources, poly_coef_matrix);
            time_this("Filling the polynomial coefficient matrix");
            double amplitude_change = Core::NaN,
                   max_amplitude_change = 0,
                   convergence_rate = 1;
            int iteration = 0;
            num_discarded_sources = 0;
            while(
                !(amplitude_change <= max_amplitude_change)
                &&
                !(convergence_rate < min_convergence_rate)
                &&
                (max_iterations<0 || iteration<max_iterations)
                &&
                num_discarded_sources == 0
            ) {
                double old_amplitude_change = amplitude_change;
                amplitude_change = fit_piecewise_bicubic_psf_step(
                        fit_sources,
                        overlap_groups,
                        symmetrized_pix_integral_matrix,
                        poly_coef_matrix,
                        modified_pixel_excesses,
                        modified_rhs_offset,
                        smoothing,
                        best_fit_poly_coef
                );
                double amplitude_vector_norm = 0;
                for(
                    LinearSourceList::iterator si=fit_sources.begin();
                    si != fit_sources.end();
                    ++si
                ) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Source ("
                              << *si
                              << ") at ("
                              << (*si)->x()
                              << ", "
                              << (*si)->y()
                              << ") flux = "
                              << (*si)->flux(0).value()
                              << std::endl;
#endif
                    amplitude_vector_norm += std::pow((*si)->flux(0).value(),2);
                }
                amplitude_vector_norm = std::sqrt(amplitude_vector_norm);
                max_amplitude_change = std::max(
                        max_abs_amplitude_change,
                        max_rel_amplitude_change * amplitude_vector_norm
                );
                convergence_rate = (
                    (old_amplitude_change - amplitude_change)
                    /
                    (old_amplitude_change - max_amplitude_change)
                );
#ifdef TRACK_PROGRESS
                std::cerr << "Amplitude vector norm="
                          << amplitude_vector_norm
                          << std::endl;
                std::cerr << "Old amplitude change=" << old_amplitude_change
                          << std::endl;
                std::cerr << "Maximum amplitude change="
                          << max_amplitude_change
                          << std::endl;
#endif
                ++iteration;
                num_discarded_sources = discard_sources(fit_sources,
                                                        max_chi2,
                                                        false,
                                                        dropped_sources);
            }
#ifdef TRACK_PROGRESS
            std::cerr << "Broke loop with iteration = " << iteration
                      << ", amplitude_change = " << amplitude_change
                      << ", max_amplitude_change = " << max_amplitude_change
                      << ", convergence_rate = " << convergence_rate
                      << ", num_discarded_sources = "
                      << num_discarded_sources << std::endl;
#endif
            if(
                convergence_rate < min_convergence_rate
                ||
                (max_iterations > 0 && iteration >= max_iterations)
            ) {
#ifdef TRACK_PROGRESS
                std::cerr << "Non convergence: convergence_rate="
                          << convergence_rate << ", min_convergence_rate="
                          << min_convergence_rate << ", iteration="
                          << iteration << ", max_iterations=" << max_iterations
                          << std::endl;
#endif
                converged = false;
                break;
            }
            if(num_discarded_sources == 0) {
                discarded_pixels = discard_outlier_pixels(
                    fit_sources,
                    pixel_integral_matrix,
                    pixel_excesses,
                    rhs_offset,
                    max_residual,
                    best_fit_poly_coef
                );
                discard_by_chi2 = true;
                num_discarded_sources += discard_sources(fit_sources,
                                                         Core::Inf,
                                                         discard_by_chi2,
                                                         dropped_sources);
            } else discarded_pixels = 0;
        }
        if(!dropped_sources.empty())
            fit_dropped_sources(dropped_sources,
                                best_fit_poly_coef);

        transform_best_fit_to_direct_parameters(basis, best_fit_poly_coef);

        convert_flux_to_adu( fit_sources, gain );
        convert_flux_to_adu( dropped_sources, gain );
        return converged;
    }

    bool fit_piecewise_bicubic_psf(LinearSourceList &fit_sources,
                                   LinearSourceList &dropped_sources,
                                   double gain,
                                   const std::vector<double> &x_grid,
                                   const std::vector<double> &y_grid,
                                   const Core::SubPixelMap &subpix_map,
                                   double max_abs_amplitude_change,
                                   double max_rel_amplitude_change,
                                   double max_chi2,
                                   double pixel_rejection,
                                   double min_convergence_rate,
                                   int max_iterations,
                                   double smoothing_penalty,
                                   Eigen::VectorXd &best_fit_poly_coef)
    {
        time_this();
        OverlapGroupList overlap_groups;
        find_overlap_groups(fit_sources, dropped_sources, overlap_groups);
        time_this("Find shape-fit overlap groups");
        estimate_initial_amplitudes(fit_sources, gain);
        time_this("Initial flux estimate");

        return fit_piecewise_bicubic_psf(fit_sources,
                                         overlap_groups,
                                         dropped_sources,
                                         x_grid,
                                         y_grid,
                                         subpix_map,
                                         max_abs_amplitude_change,
                                         max_rel_amplitude_change,
                                         max_chi2,
                                         pixel_rejection,
                                         min_convergence_rate,
                                         max_iterations,
                                         smoothing_penalty,
                                         best_fit_poly_coef,
                                         gain);
    }

    bool fit_piecewise_bicubic_psf(LinearSourceList &fit_sources,
                                   LinearSourceList &dropped_sources,
                                   const PSF::PiecewiseBicubicMap &psf_map,
                                   const std::vector<double> &x_grid,
                                   const std::vector<double> &y_grid,
                                   const Core::SubPixelMap &subpix_map,
                                   double max_abs_amplitude_change,
                                   double max_rel_amplitude_change,
                                   double max_chi2,
                                   double pixel_rejection,
                                   double min_convergence_rate,
                                   int max_iterations,
                                   double smoothing_penalty,
                                   Eigen::VectorXd &best_fit_poly_coef,
                                   double gain)
    {
        time_this();
        OverlapGroupList overlap_groups;
        find_overlap_groups(fit_sources, dropped_sources, overlap_groups);
        time_this("Find shape-fit overlap groups");
        estimate_initial_amplitudes(fit_sources,
                                    overlap_groups,
                                    psf_map);
        time_this("Initial flux estimate");

        unsigned expansion_rows = ((psf_map.x_grid().size() - 2)
                                   *
                                   (psf_map.y_grid().size() - 2)),
                 num_coef = 4 * expansion_rows * psf_map.num_terms();
        best_fit_poly_coef.resize(num_coef);
        best_fit_poly_coef = PSF::RowMajorMap(psf_map.coefficients(),
                                              num_coef,
                                              1);
        if(max_iterations == 0) {
            if(!dropped_sources.empty())
                fit_dropped_sources(dropped_sources,
                                    best_fit_poly_coef);
            return true;
        } else
            return fit_piecewise_bicubic_psf(fit_sources,
                                             overlap_groups,
                                             dropped_sources,
                                             x_grid,
                                             y_grid,
                                             subpix_map,
                                             max_abs_amplitude_change,
                                             max_rel_amplitude_change,
                                             max_chi2,
                                             pixel_rejection,
                                             min_convergence_rate,
                                             max_iterations,
                                             smoothing_penalty,
                                             best_fit_poly_coef,
                                             gain);
    }

    void fit_dropped_sources(
        LinearSourceList &dropped_sources,
        const Eigen::VectorXd &best_fit_poly_coef
    )
    {
        for(
            LinearSourceList::iterator si = dropped_sources.begin();
            si != dropped_sources.end();
            ++si
        ) {
            if((*si)->overlaps().size() == 0) {
                (*si)->fit_flux(best_fit_poly_coef);
            }
#ifdef VERBOSE_DEBUG
            std::cerr << "Source ("
                      << *si
                      << ") at ("
                      << (*si)->x()
                      << ", "
                      << (*si)->y()
                      << ") flux = "
                      << (*si)->flux(0).value()
                      << std::endl;
#endif
        }
    }

    void output_best_fit_psf(const Eigen::VectorXd &best_fit_poly_coef,
                             const std::vector<double> &x_grid,
                             const std::vector<double> &y_grid,
                             const std::string &fname)
    {
#ifdef DEBUG
        assert(best_fit_poly_coef.size() % (x_grid.size() * y_grid.size()) == 0);
#endif
        std::ofstream outf(fname.c_str());
        outf.setf(std::ios_base::scientific);
        outf.precision(16);
        outf << "x grid " << std::endl << x_grid.size() << ":";
        for(Core::vector_size_type i = 0; i < x_grid.size(); ++i)
            outf << " " << x_grid[i];
        outf << std::endl << "y grid " << std::endl << y_grid.size() << ":";
        for(Core::vector_size_type i = 0; i < y_grid.size(); ++i)
            outf << " " << y_grid[i];
        outf << std::endl << "coefficients:";
        for(int i = 0; i < best_fit_poly_coef.size(); ++i)
            outf << " " << best_fit_poly_coef[i];
        outf << std::endl;
        outf.close();
    }

} //End FitPSF namespace.
