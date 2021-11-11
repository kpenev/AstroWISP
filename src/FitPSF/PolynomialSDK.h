/**\file
 *
 * \brief Declares a class for fitting elliptical gaussian PSFs as a smooth
 * function of image position
 *
 * \ingroup FitPSF
 */

#ifndef __POLYNOMIAL_SDK_H
#define __POLYNOMIAL_SDK_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SDKUtil.h"
#include "../PSF/EllipticalGaussian.h"
#include "../PSF/MapSource.h"
#include "../IO/OutputSDKSource.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <cmath> 

namespace FitPSF {

    ///\brief A class that assumes the PSF on an image is described by an
    ///elliptical Gaussian, specified via the S, D, K parameters, and that
    ///those parameters are low order polynomials of the image position and
    ///fits for the coefficients of the polynomials.
    template< class SUBPIX_TYPE >
        class LIB_LOCAL PolynomialSDK {
        private:
            double
                ///Lower limit for the initial search of PSF scale
                __minS,

                ///Upper limit for the initial search of PSF scale
                __maxS,

                ///The tolerance to which to determine fit values
                __fit_tolerance,

                ///The maximum reduced chi squared allowed before a source is
                ///dropped from the fit.
                __max_source_chi2;

            ///\brief The polynomial coefficients of S (first), D (second), K
            ///(third) in the order of the PSFMap terms.
            std::valarray<double> __poly_coef;


            ///The largest absolute value of each PSF term over all sources.
            Eigen::VectorXd __max_abs_expansion_terms;

            ///\brief Returns initial guesses for the S, D and K polynomials
            ///in a newly allocated GSL vector generated from the
            ///coefficients of a fit with one less term.
            gsl_vector *initial_poly(gsl_vector *lower_order_coef);

            ///Finds a point between min_S and max_S that has a chi squared
            ///value less than at both min_S and max_S. Also possibly updates
            ///min_S and max_S to a narrower internal.
            double find_mid_S(double &min_S,
                              double &max_S,
                              void *param_array);

            ///Performs a sigle dimensional fit for an approximate S value of
            ///the image.
            template<class SOURCE_ITERATOR>
                double fit_initial_S(SOURCE_ITERATOR first_source,
                                     SOURCE_ITERATOR past_last_source,
                                     const SUBPIX_TYPE &subpix_map,
                                     double min_S,
                                     double max_S);

            ///Performs a fit for the PSF dependence restricted to terms of
            ///up to the given order starting with the given initial
            ///coefficients. Assumes that gsl_params variable is already
            ///properly initialized as required by gsl_minimization_function.
            void restricted_gsl_fit(
                gsl_vector* poly_coef,
                void*       gsl_params,
                double      initial_step_size = 0.5
            );

            ///Performs a single nr step (updating poly_coef and chi2) such
            ///that at the end chi squared is smaller than at the beginning.
            template<class SOURCE_ITERATOR>
                bool nr_step(
                    SOURCE_ITERATOR first_source,
                    SOURCE_ITERATOR past_last_source,
                    const SUBPIX_TYPE &subpix_map,
                    const Eigen::VectorXd &coef_change,
                    double max_SDK_change,
                    Eigen::VectorXd &poly_coef, double &chi2,
                    Eigen::VectorXd &grad, Eigen::MatrixXd &d2
                );

            ///Performs a Newton-Raphson fit for a single order.
            template<class SOURCE_ITERATOR>
                void single_order_nr_fit(
                    SOURCE_ITERATOR first_source,
                    SOURCE_ITERATOR past_last_source,
                    const SUBPIX_TYPE &subpix_map,
                    Eigen::VectorXd &poly_coef
                );

            ///Fills in __max_abs_expansion_terms.
            template<class SOURCE_ITERATOR>
                void find_max_expansion_terms(
                    SOURCE_ITERATOR first_source,
                    SOURCE_ITERATOR past_last_source
                );
        public:
            ///Construct a polynomial S, D, K dependence for fitting.
            PolynomialSDK(
                ///The lower end of the range of S values to consider for
                ///this image.
                double minS,

                ///The upper end of the range of S values to consider for
                ///this image.
                double maxS,

                ///The maximum allowed uncertainty in the best fit values of
                ///S, D and K
                double fit_tolerance = 1e-2,

                ///Sources for which the reduced chi2 is larger than this
                ///value are dropped from the fit on the grounds that they
                ///are probably not point sources.
                double max_source_chi2 = 10.0
            ) :
                __minS(minS),
                __maxS(maxS),
                __fit_tolerance(fit_tolerance),
                __max_source_chi2(max_source_chi2)
            {}

            ///\brief Performs the fit for the PSF dependence based on the
            ///given list of sources and returns a list of SDK sources with
            ///properly set S,D,K, amplitude and background. Uses GSL simplex
            ///method for the fit.
            template<class SOURCE_ITERATOR>
                void gsl_fit(SOURCE_ITERATOR first_source,
                             SOURCE_ITERATOR past_last_source,
                             const SUBPIX_TYPE &subpix_map,
                             const std::list<double>
                             &guess_poly_coef = std::list<double>());

            ///\brief Fills grad with the gradient vector and d2 with the
            ///second derivative matrix of chi squared for the given sources
            ///and S, D, K expansion of the given order. The gradient (second
            ///order derivative) is only defined if the given sources have
            ///first (second) order derivatives enabled.
            template<class SOURCE_ITERATOR>
                void chi2_grad_d2(
                    SOURCE_ITERATOR first_source,
                    SOURCE_ITERATOR past_last_source,
                    const Eigen::VectorXd &poly_coef,
                    const SUBPIX_TYPE &subpix_map, double &chi2,
                    Eigen::VectorXd &grad, Eigen::MatrixXd &d2
                );

            ///\brief The same as gsl_fit, but uses Newton-Raphon method to
            ///solve for a zero of the gradient of chi squared.
            template<class SOURCE_ITERATOR>
                void nr_fit(SOURCE_ITERATOR first_source,
                            SOURCE_ITERATOR past_last_source,
                            const SUBPIX_TYPE &subpix_map,
                            const std::list<double>
                            &guess_poly_coef=std::list<double>());

            ///\brief Prepares a list of SDKSource objects properly
            ///initialized with S, D, K, amplitude and background according
            ///to the latest best fit polynomial coefficients.
            template<class SOURCE_ITERATOR>
                std::list<IO::OutputSDKSource> best_fit_sources(
                    SOURCE_ITERATOR first_source,
                    SOURCE_ITERATOR past_last_source,
                    const SUBPIX_TYPE &subpix_map,
                    double gain
                );

            ///Make PolynomialSDK::operator() available.
            template<class COEF_TYPE>
                PSF::EllipticalGaussian operator()(
                    const PSF::MapSource &source,
                    const COEF_TYPE &poly_coef
                );

            ///\brief Evaluates the most recently derived fit (or the
            ///polynomial defined by the given coefficients) at the given
            ///position, returning a PSF that would be appropriate for a
            ///source centered there.
            PSF::EllipticalGaussian operator()(const PSF::MapSource &source)
            {return operator()(source, __poly_coef);}

            ///Returns the best fit polynomial coefficients.
            const std::valarray<double> &get_coefficients() const
            {return __poly_coef;}
        }; //End PolynomialSDK class.

    ///\brief Returns initial guesses for the S, D and K polynomials in a
    ///newly allocated GSL vector generated from the coefficients of a fit
    ///one smaller order.
    template<class SUBPIX_TYPE>
        gsl_vector *PolynomialSDK<SUBPIX_TYPE>::initial_poly(
            gsl_vector *lower_order_coef
        )
        {
            if(lower_order_coef->size%3) throw Error::InvalidArgument(
                "PolynomialSDK::initial_poly",
                "Previous polynomial coefficient count not divisible by 3."
            );
            unsigned known_coef_count=lower_order_coef->size/3,
                     new_coef_count=known_coef_count + 1;
            gsl_vector *result=gsl_vector_calloc(3*new_coef_count);

            for(unsigned var_ind=0; var_ind<3; var_ind++)
                for(unsigned i=0; i<known_coef_count; i++)
                    gsl_vector_set(
                        result,
                        i + var_ind * new_coef_count,
                        gsl_vector_get(
                            lower_order_coef,
                            i + var_ind *known_coef_count
                        )
                    );
            return result;
        }

    ///\brief Finds a point between min_S and max_S that has a chi squared
    ///value less than at both min_S and max_S. Also possibly updates min_S
    ///and max_S to a narrower internal.
    template<class SUBPIX_TYPE>
        double PolynomialSDK<SUBPIX_TYPE>::find_mid_S(double &min_S,
                                                      double &max_S,
                                                      void *gsl_param)
        {
            std::list<double> test_s,
                              chi2_values;
            double chi2_mid, s_mid;
            test_s.push_back(min_S);
            chi2_values.push_back(gsl_s_minimization_function(min_S,
                                                              gsl_param));
            test_s.push_back(max_S);
            chi2_values.push_back(gsl_s_minimization_function(max_S,
                                                              gsl_param));
            std::list<double>::iterator s_right = test_s.begin(),
                                        s_left = s_right++,
                                        chi2_right = chi2_values.begin(),
                                        chi2_left = chi2_right++;
            bool searching;
            do {
                if(*s_right - *s_left < 10.0 * __fit_tolerance) {
                    if(chi2_values.front() < chi2_values.back()) {
                        max_S = min_S;
                        return min_S;
                    } else {
                        min_S = max_S;
                        return max_S;
                    }
                }
                s_mid = 0.5 * (*s_left + *s_right);
                chi2_mid = gsl_s_minimization_function(s_mid, gsl_param);
                searching = (chi2_mid >= *chi2_left
                             ||
                             chi2_mid >= *chi2_right);
                if(searching) {
                    test_s.insert(s_right, s_mid);
                    chi2_values.insert(chi2_right, chi2_mid);
                    s_left = s_right++;
                    if(s_right == test_s.end()) {
                        s_right = test_s.begin();
                        s_left = s_right++;
                        chi2_right = chi2_values.begin();
                    }
                    chi2_left = chi2_right++;
                }
            } while (searching);
            min_S = *s_left;
            max_S = *s_right;
            return s_mid;
        }

    ///Performs a sigle dimensional fit for an approximate S value of the
    ///image.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        double PolynomialSDK<SUBPIX_TYPE>::fit_initial_S(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            double min_S,
            double max_S
        )
    {
        void *param_array[3];
        param_array[0] = reinterpret_cast<void*>(&first_source);
        param_array[1] = reinterpret_cast<void*>(&past_last_source);
        param_array[2] = reinterpret_cast<void*>(
            const_cast<SUBPIX_TYPE*>(&subpix_map)
        );
        void *gsl_param = reinterpret_cast<void*>(param_array);

        const gsl_min_fminimizer_type *minimizer_type =
            gsl_min_fminimizer_brent;
        double mid_S = find_mid_S(min_S, max_S, gsl_param);
        if(max_S - min_S < __fit_tolerance) return mid_S;

        gsl_function minimization_function;

        minimization_function.function = &gsl_s_minimization_function;
        minimization_function.params = gsl_param;
        gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(
            minimizer_type
        );
        gsl_min_fminimizer_set(minimizer,
                               &minimization_function,
                               mid_S,
                               min_S,
                               max_S);

        int iter = 0;

        for(
            int gsl_status = GSL_CONTINUE;
            gsl_status == GSL_CONTINUE; ++iter
        ) {
            gsl_status = gsl_min_fminimizer_iterate(minimizer);

            mid_S = gsl_min_fminimizer_x_minimum(minimizer);
            min_S = gsl_min_fminimizer_x_lower(minimizer);
            max_S = gsl_min_fminimizer_x_upper(minimizer);

            gsl_status = gsl_min_test_interval(min_S, max_S, 0.01, 0.0);
        }

        gsl_min_fminimizer_free(minimizer);
        return mid_S;
    }

    ///\brief Performs a fit for the PSF dependence restricted to terms of up
    ///to the given order starting with the given initial coefficients
    ///(overwritten on exit). Assumes that gsl_params variable is already
    ///properly initialized as required by gsl_minimization_function.
    template<class SUBPIX_TYPE>
        void PolynomialSDK<SUBPIX_TYPE>::restricted_gsl_fit(
            gsl_vector* poly_coef,
            void*       gsl_params,
            double      initial_step_size
        )
    {
        const gsl_multimin_fminimizer_type *MinimizerType =
            gsl_multimin_fminimizer_nmsimplex;
        gsl_multimin_fminimizer *minimizer;
        gsl_vector *step_size = gsl_vector_alloc(poly_coef->size);
        gsl_vector_set_all(step_size, initial_step_size);
        gsl_multimin_function minimizer_func;
        minimizer_func.n = poly_coef->size;
        minimizer_func.f = gsl_minimization_function;
        minimizer_func.params = gsl_params;
        minimizer = gsl_multimin_fminimizer_alloc(MinimizerType,
                                                  poly_coef->size);
        gsl_multimin_fminimizer_set(minimizer, &minimizer_func, poly_coef,
                step_size);

        int status = GSL_CONTINUE;
        for(unsigned iter = 0; status == GSL_CONTINUE; ++iter) {
            status = gsl_multimin_fminimizer_iterate(minimizer);

            if (status)
                throw Error::GSLError("Nonzero return status from iterate "
                        "while minimizing PSF map.");

            double size = gsl_multimin_fminimizer_size (minimizer);
            status = gsl_multimin_test_size (size, __fit_tolerance);

            if (status == GSL_SUCCESS) break;
        }
        for(unsigned i = 0; i < poly_coef->size; ++i)
            gsl_vector_set(poly_coef, i, gsl_vector_get(minimizer->x, i));
        gsl_vector_free(step_size);
        gsl_multimin_fminimizer_free (minimizer);
    }

    ///Performs a single nr step (updating poly_coef and chi2) such that at
    ///the end chi squared is smaller than at the beginning.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        bool PolynomialSDK<SUBPIX_TYPE>::nr_step(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            const Eigen::VectorXd &coef_change,
            double max_SDK_change,
            Eigen::VectorXd &poly_coef,
            double &chi2,
            Eigen::VectorXd &grad,
            Eigen::MatrixXd &d2)
        {
#ifdef DEBUG
            assert(!std::isnan(coef_change.sum()));
#endif
            double gamma=1, orig_chi2=chi2;
            Eigen::VectorXd new_poly_coef;
            while(orig_chi2 <= chi2 || std::isnan(chi2)) {
                double new_s = poly_coef[0] + gamma * coef_change[0];
                if(poly_coef.size() <= 3) {
                    if(new_s < __minS)
                        gamma = (__minS - poly_coef[0]) / coef_change[0];
                    else if(new_s>__maxS)
                        gamma = (__maxS - poly_coef[0]) / coef_change[0];
                    if(poly_coef.size() == 3) {
                        double s2_m_d2_m_k2 = (
                            std::pow(poly_coef[0], 2)
                            -
                            std::pow(poly_coef[1], 2)
                            -
                            std::pow(poly_coef[2], 2)
                        );
                        if(s2_m_d2_m_k2 < 0) {
                            double
                                a = (std::pow(coef_change[0], 2)
                                     -
                                     std::pow(coef_change[1], 2)
                                     -
                                     std::pow(coef_change[2], 2)),
                                b = 2.0 * (poly_coef[0] * coef_change[0]
                                           -
                                           poly_coef[1] * coef_change[1]
                                           -
                                           poly_coef[2] * coef_change[2]),
                                sol1,
                                sol2;
                            int nroots = gsl_poly_solve_quadratic(
                                a,
                                b, s2_m_d2_m_k2,
                                &sol1,
                                &sol2
                            );
                            if(nroots == 1) {
#ifdef DEBUG
                                assert(sol1 * gamma>0);
#endif
                                gamma = sol1;
                            } else if(nroots == 2) {
                                if(sol1 * gamma > 0) gamma = sol1;
                                else if(sol2 * gamma > 0) gamma = sol2;
#ifdef DEBUG
                                else assert(false);
#endif
                            }
                        }
                    }
                }
#ifdef DEBUG
                assert(!std::isnan(poly_coef.sum()));
                assert(!std::isnan(gamma));
#endif
                new_poly_coef = poly_coef + gamma * coef_change;
                chi2_grad_d2(first_source,
                             past_last_source,
                             new_poly_coef,
                             subpix_map,
                             chi2,
                             grad,
                             d2);
                if(std::abs(gamma) > __fit_tolerance / max_SDK_change)
                    gamma /= 10.0;
                else if(gamma > 0) gamma = -1;
                else {
                    chi2 = orig_chi2;
                    return false;
                }
            }
            poly_coef = new_poly_coef;
            return true;
        }

    ///Performs a Newton-Raphson fit for a single order.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        void PolynomialSDK<SUBPIX_TYPE>::single_order_nr_fit(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            Eigen::VectorXd &poly_coef
        )
        {
            double chi2,
                   max_SDK_change;
            Eigen::VectorXd grad,
                            coef_change,
                            max_poly_terms(poly_coef.size());
            Eigen::MatrixXd d2;

            chi2_grad_d2(first_source,
                         past_last_source,
                         poly_coef,
                         subpix_map,
                         chi2,
                         grad,
                         d2);
#ifdef DEBUG
            assert(!std::isnan(grad.sum()));
            assert(!std::isnan(d2.sum()));
            assert(poly_coef.size() == 1 || poly_coef.size() % 3 == 0);
#endif

            unsigned fit_nterms = std::max(long(1), poly_coef.size() / 3);

            while(true) {
                if(poly_coef.size() == 1) {
                    coef_change.resize(1);
                    coef_change[0] = -grad[0] / d2(0,0);
                } else {
                    Eigen::FullPivHouseholderQR<Eigen::MatrixXd>
                        decomposition(d2);
                    coef_change = decomposition.solve(-grad);
                }
                max_SDK_change = 0;
                for(
                    unsigned var_i = 0;
                    var_i < (poly_coef.size() == 1 ? 1 : 3);
                    ++var_i
                )
                    max_SDK_change = std::max(
                        max_SDK_change,
                        coef_change.segment(
                            var_i * fit_nterms, fit_nterms
                        ).cwiseAbs().cwiseProduct(
                            __max_abs_expansion_terms.head(fit_nterms)
                        ).maxCoeff()
                    );
#ifdef TRACK_PROGRESS
                std::cerr << "Maximum SDK change = " << max_SDK_change
                    << " (fit tolerance = " << __fit_tolerance << ")."
                    << std::endl;
#endif
                if(max_SDK_change < __fit_tolerance) break;
                if(
                    ! nr_step(first_source,
                              past_last_source,
                              subpix_map,
                              coef_change,
                              max_SDK_change,
                              poly_coef,
                              chi2,
                              grad,
                              d2)
                ) {
                    double grad_scale = grad.cwiseAbs().cwiseProduct(
                        __max_abs_expansion_terms.head(grad.size())
                    ).maxCoeff() / max_SDK_change;
                    coef_change = -grad / grad_scale;
                    if(
                        !nr_step(first_source,
                                 past_last_source,
                                 subpix_map,
                                 coef_change,
                                 max_SDK_change,
                                 poly_coef,
                                 chi2,
                                 grad,
                                 d2)
                    ) return;
                }
            }

        }

    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        void PolynomialSDK<SUBPIX_TYPE>::find_max_expansion_terms(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source
        )
        {
            unsigned nterms = first_source->expansion_terms().size();
            __max_abs_expansion_terms = Eigen::VectorXd::Zero(nterms);
            for(
                SOURCE_ITERATOR si = first_source;
                si != past_last_source;
                ++si
            ) {
                const Eigen::VectorXd
                    &expansion_terms = si->expansion_terms();
                for(unsigned term_i = 0; term_i < nterms; ++term_i)
                    __max_abs_expansion_terms[term_i] = std::max(
                        std::abs(expansion_terms[term_i]),
                        __max_abs_expansion_terms[term_i]
                    );
            }
        }

    ///\brief Fills grad with the gradient vector and d2 with the second
    ///derivative matrix of chi squared for the given sources and S, D, K
    ///expansion of the given order. The gradient (second order derivative)
    ///is only defined if the given sources have first (second) order
    ///derivatives enabled.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        void PolynomialSDK<SUBPIX_TYPE>::chi2_grad_d2(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const Eigen::VectorXd &poly_coef,
            const SUBPIX_TYPE &subpix_map,
            double &chi2,
            Eigen::VectorXd &grad,
            Eigen::MatrixXd &d2
        )
        {
            unsigned coef_count = poly_coef.size() / 3;

            chi2 = 0.0;
            if(coef_count) {
                grad = Eigen::VectorXd::Zero(3 * coef_count);
                d2 = Eigen::MatrixXd::Zero(3 * coef_count, 3 * coef_count);
            } else {
                grad = Eigen::VectorXd::Zero(1);
                d2 = Eigen::MatrixXd::Zero(1, 1);
            }

            for(
                SOURCE_ITERATOR si = first_source;
                si != past_last_source;
                ++si
            ) {
                if(si->is_nonpoint()) continue;
                PSF::EllipticalGaussian psf = operator()(*si, poly_coef);
                if(
                    psf.s() < __minS || psf.s() > __maxS
                    ||
                    std::pow(psf.s(), 2) <= std::pow(psf.d(), 2)
                    +
                    std::pow(psf.k(), 2)
                ) {
                    chi2 = Inf;
                    return;
                }

                assert(!std::isnan(psf.s()));
                assert(!std::isnan(psf.d()));
                assert(!std::isnan(psf.k()));

                si->set_PSF(psf, subpix_map);
                const std::valarray<double>&
                    chi2_deriv = si->chi2_all_deriv();

                assert(!std::isnan(chi2_deriv[PSF::NO_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::S_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::D_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::K_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::SS_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::SD_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::SK_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::DD_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::DK_DERIV]));
                assert(!std::isnan(chi2_deriv[PSF::KK_DERIV]));

                chi2 += chi2_deriv[PSF::NO_DERIV];

                if(coef_count) {
                    grad.head(coef_count) += (
                        si->expansion_terms().head(coef_count)
                        *
                        chi2_deriv[PSF::S_DERIV]
                    );
                    grad.segment(coef_count, coef_count) += (
                        si->expansion_terms().head(coef_count)
                        *
                        chi2_deriv[PSF::D_DERIV]
                    );
                    grad.tail(coef_count) += (
                        si->expansion_terms().head(coef_count)
                        *
                        chi2_deriv[PSF::K_DERIV]
                    );

                    Eigen::MatrixXd two_poly_terms = (
                        si->expansion_terms().head(coef_count)
                        *
                        si->expansion_terms().head(coef_count).transpose()
                    );

                    d2.topLeftCorner(coef_count, coef_count) +=
                        two_poly_terms * chi2_deriv[PSF::SS_DERIV];
                    d2.block(0, coef_count, coef_count, coef_count) +=
                        two_poly_terms * chi2_deriv[PSF::SD_DERIV];
                    d2.block(0, 2*coef_count, coef_count, coef_count) +=
                        two_poly_terms * chi2_deriv[PSF::SK_DERIV];
                    d2.block(coef_count, coef_count, coef_count, coef_count) +=
                        two_poly_terms * chi2_deriv[PSF::DD_DERIV];
                    d2.block(coef_count,
                             2*coef_count,
                             coef_count,
                             coef_count) += (two_poly_terms
                                             *
                                             chi2_deriv[PSF::DK_DERIV]);
                    d2.bottomRightCorner(coef_count, coef_count) +=
                        two_poly_terms * chi2_deriv[PSF::KK_DERIV];
                } else {
                    grad(0) += chi2_deriv[PSF::S_DERIV];
                    d2(0,0) += chi2_deriv[PSF::SS_DERIV];
                }
            }
            if(coef_count) {
                d2.block(coef_count, 0, coef_count, coef_count) =
                    d2.block(0, coef_count, coef_count, coef_count);

                d2.block(2*coef_count, 0, coef_count, coef_count) =
                    d2.block(0, 2*coef_count, coef_count, coef_count);

                d2.block(2*coef_count, coef_count, coef_count, coef_count) =
                    d2.block(coef_count,
                             2*coef_count,
                             coef_count,
                             coef_count);
            }
        }

    ///Performs the fit for the PSF dependence based on the given list of
    ///sources and returns a list of SDK sources with properly set S,D,K,
    ///amplitude and background. Uses GSL simplex method for the fit.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        void PolynomialSDK<SUBPIX_TYPE>::gsl_fit(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            const std::list<double> &guess_poly_coef
        )
        {
            find_max_expansion_terms(first_source, past_last_source);
            void *param_array[7];
            double max_source_chi2=Inf;
            param_array[0] = reinterpret_cast<void*>(&first_source);
            param_array[1] = reinterpret_cast<void*>(&past_last_source);
            param_array[2] = reinterpret_cast<void*>(this);
            param_array[3] = reinterpret_cast<void*>(
                const_cast<SUBPIX_TYPE*>(&subpix_map));
            param_array[4] = reinterpret_cast<void*>(&max_source_chi2);
            param_array[5] = reinterpret_cast<void*>(&__minS);
            param_array[6] = reinterpret_cast<void*>(&__maxS);
            void *gsl_params = reinterpret_cast<void*>(param_array);
            gsl_vector *poly_coef;
            if(guess_poly_coef.size() == 0) {
                poly_coef = gsl_vector_calloc(3);
                gsl_vector_set(poly_coef, 0, fit_initial_S(first_source,
                                                           past_last_source,
                                                           subpix_map,
                                                           __minS, __maxS));
                for(
                    unsigned fit_nterms = 0;
                    fit_nterms <= __max_abs_expansion_terms.size();
                    ++fit_nterms
                ) {
                    restricted_gsl_fit(poly_coef,
                                       gsl_params,
                                       10.0 * __fit_tolerance );
                    if(fit_nterms < __max_abs_expansion_terms.size()) {
                        gsl_vector *new_poly_coef = initial_poly(poly_coef);
                        gsl_vector_free(poly_coef);
                        poly_coef = new_poly_coef;
                    }
                }
            } else {
#ifdef DEBUG
                assert(static_cast<int>(guess_poly_coef.size())
                       ==
                       static_cast<int>(__max_abs_expansion_terms.size()));
#endif
                poly_coef = gsl_vector_alloc(guess_poly_coef.size());
                std::list<double>::const_iterator
                    coef_iter = guess_poly_coef.begin();
                for(unsigned i = 0; i < guess_poly_coef.size(); i++)
                    gsl_vector_set(poly_coef, i, *coef_iter++);
            }
            max_source_chi2 = __max_source_chi2;
            restricted_gsl_fit(poly_coef, gsl_params, 2.0 * __fit_tolerance);
            __poly_coef.resize(poly_coef->size);
            for(unsigned i = 0; i < poly_coef->size; i++)
                __poly_coef[i] = gsl_vector_get(poly_coef, i);
        }

    ///The same as gsl_fit, but uses Newton-Raphon method to solve for a zero
    ///of the gradient of chi squared.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        void PolynomialSDK<SUBPIX_TYPE>::nr_fit(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            const std::list<double> &guess_poly_coef
        )
        {
            find_max_expansion_terms(first_source, past_last_source);
            Eigen::VectorXd poly_coef;

            if(guess_poly_coef.size() == 0) {
                poly_coef.setConstant(1, 0.5 * (__minS + __maxS));
            } else {
#ifdef DEBUG
                assert(static_cast<int>(guess_poly_coef.size())
                       ==
                       3 * static_cast<int>(__max_abs_expanson_terms.size()));
#endif
                poly_coef.resize(guess_poly_coef.size());
                std::list<double>::const_iterator
                    coef_iter = guess_poly_coef.begin();
                for(unsigned i = 0; i < guess_poly_coef.size(); i++)
                    poly_coef[i] = *coef_iter++;
            }

            bool done = false;//guess_poly_coef.size() > 0;
            while(!done) {
#ifdef TRACK_PROGRESS
                std::cerr << "Performing " << poly_coef.size()
                    << " coefficient fit." << std::endl;
#endif
                single_order_nr_fit(first_source,
                                    past_last_source,
                                    subpix_map,
                                    poly_coef);
                done = poly_coef.size() == (
                    3 * __max_abs_expansion_terms.size()
                );
                if(!done) {
                    Eigen::VectorXd old_poly_coef=poly_coef;
                    unsigned old_coef_count = poly_coef.size() / 3;
                    unsigned new_coef_count = old_coef_count + 1;
                    poly_coef.setZero(3 * new_coef_count);
                    if(old_coef_count == 0) poly_coef[0] = old_poly_coef[0];
                    else {
                        poly_coef.head(old_coef_count)=
                            old_poly_coef.head(old_coef_count);

                        poly_coef.segment(new_coef_count, old_coef_count)=
                            old_poly_coef.segment(old_coef_count,
                                                  old_coef_count);

                        poly_coef.segment(2 * new_coef_count,
                                          old_coef_count) =
                            old_poly_coef.tail(old_coef_count);
                    }
                }
            }

        unsigned num_rejected = 1;
        while(num_rejected > 0) {
            num_rejected = 0;
            for(
                SOURCE_ITERATOR si = first_source;
                si != past_last_source;
                ++si
            )
                if(si->chi2() / (si->fit_pixels() - 1) > __max_source_chi2)
                    if(!si->is_nonpoint()) {
                        si->set_nonpoint();
                        ++num_rejected;
                    }
            if(num_rejected > 0)
                single_order_nr_fit(
                        first_source, past_last_source, subpix_map, poly_coef
                );
        }

        __poly_coef.resize(poly_coef.size());
        for(int i = 0; i < poly_coef.size(); i++)
            __poly_coef[i] = poly_coef[i];
    }

    ///Prepares a list of SDKSource objects properly initialized with S, D,
    ///K, amplitude and background according to the latest best fit
    ///polynomial coefficients.
    template<class SUBPIX_TYPE>
        template<class SOURCE_ITERATOR>
        std::list<IO::OutputSDKSource>
        PolynomialSDK<SUBPIX_TYPE>::best_fit_sources(
            SOURCE_ITERATOR first_source,
            SOURCE_ITERATOR past_last_source,
            const SUBPIX_TYPE &subpix_map,
            double gain)
        {
            std::list<IO::OutputSDKSource> result;
            int source_id=1;
            for(GSLSourceIteratorType src_iter=first_source;
                src_iter!=past_last_source; src_iter++) {
                PSF::EllipticalGaussian psf=(*this)(*src_iter);
                src_iter->set_PSF(psf, subpix_map);
                double s=psf.s(), d=psf.d(), k=psf.k();
                result.push_back(
                    IO::OutputSDKSource(
                        src_iter->image_filename(),
                        src_iter->output_filename(),
                        src_iter->id(),
                        1,
                        src_iter->x(),
                        src_iter->y(),
                        s,
                        d,
                        k,
                        src_iter->amplitude() / gain,
                        SourceBackground(
                            src_iter->background_electrons() / gain,
                            std::sqrt(
                                src_iter->background_electrons_variance()
                            )
                            /
                            gain,
                            src_iter->background_pixels()
                        ),
                        src_iter->chi2() / (src_iter->fit_pixels() - 1),
                        src_iter->fit_pixels(),
                        std::sqrt(src_iter->merit())
                    )
                );

                Flux& flux = result.back().flux(0); // shortcut
                flux.value() = (
                    2 * M_PI * src_iter->amplitude()
                    /
                    std::sqrt( s*s - d*d - k*k )
                    /
                    gain
                );
                flux.error() =
                    std::sqrt(
                        flux.value() * gain
                        +
                        (
                            src_iter->background_electrons()
                            +
                            src_iter->background_electrons_variance()
                        )
                        *
                        src_iter->pixel_count()
                    ) / gain;

                flux.flag() = src_iter->quality_flag();

                ++source_id;
            }
            return result;
        }

    template<class SUBPIX_TYPE>
        template<class COEF_TYPE>
        PSF::EllipticalGaussian PolynomialSDK<SUBPIX_TYPE>::operator()(
            const PSF::MapSource &source,
            const COEF_TYPE &poly_coef
        )
        {
            double s = 0, d = 0, k = 0;
            if(poly_coef.size() == 1) s = poly_coef[0];
            else {
                assert(poly_coef.size() % 3 == 0);
                unsigned var_size = poly_coef.size() / 3;
                const Eigen::VectorXd
                    &expansion_terms = source.expansion_terms();
                for(size_t i=0; i<var_size; i++) {
                    s += expansion_terms[i]*poly_coef[i];
                    d += expansion_terms[i]*poly_coef[var_size + i];
                    k += expansion_terms[i]*poly_coef[2 * var_size + i];
                }
            }
            return PSF::EllipticalGaussian(s, d, k, 0, 1e-5, 0.0);
        }

} // End FitPSF namespace.

#endif
