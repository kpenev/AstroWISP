/**\file
 *
 * \brief The implementation of the methods of
 * EllipticalGaussianIntegralWedge.
 *
 * \ingroup PSF
 */

#include "EllipticalGaussianIntegralWedge.h"

namespace PSF {

    ///\brief Use one order higher in the given index estimates of the
    ///integral and any derivatives.
    void EllipticalGaussianIntegralWedge::update_estimates(unsigned index)
    {
        __multipliers[index].push_back(-__multipliers[index].back()
                                       *
                                       __wedge_coef[index]
                                       /
                                       __orders[index]);

        std::valarray<double> new_terms(0.0, KK_DERIV + 1);

        unsigned i_endindex=__orders[0];
        unsigned i_startindex=(index ? 0 : i_endindex);
        for(unsigned i = i_startindex; i <= i_endindex; ++i) {

            unsigned j_endindex = __orders[1];
            unsigned j_startindex = (index==1 ? j_endindex : 0);
            for(unsigned j = j_startindex; j <= j_endindex; ++j) {

                const double ij_term = (__multipliers[0][i]
                                        *
                                        __multipliers[1][j]);
                unsigned k_endindex = __orders[2];
                unsigned k_startindex = (index == 2 ? k_endindex : 0);
                for(unsigned k = k_startindex; k <= k_endindex; ++k) {

                    const double ijk_term = ij_term * __multipliers[2][k];
                    unsigned l_endindex = __orders[3];
                    unsigned l_startindex = (index == 3 ? l_endindex : 0);
                    for(unsigned l = l_startindex; l <= l_endindex; ++l) {

                        const double ijkl_term = (ijk_term
                                                  *
                                                  __multipliers[3][l]);
                        unsigned m_endindex = __orders[4];
                        unsigned m_startindex = (index==4 ? m_endindex : 0);
                        for(unsigned m = m_startindex; m <= m_endindex; ++m){
                            fill_new_terms(
                                (
                                    ijkl_term
                                    *
                                    __multipliers[4][m]
                                    *
                                    __integral(2 * i + j + l, 2 * k + j + m)
                                ),
                                i,
                                j,
                                k,
                                l,
                                m,
                                new_terms
                            );
                        }
                    }
                }
            }
        }

        __value[NO_DERIV] += __center_psf_r2 * new_terms[NO_DERIV];

        int d_startindex = (__first_deriv ? S_DERIV : SS_DERIV);
        int d_endindex = (__second_deriv ? KK_DERIV : K_DERIV);
        for(int d = d_startindex; d <= d_endindex; ++d) 
            __value[d] += __center_psf_r2 * new_terms[d];
    }

    EllipticalGaussianIntegralWedge::EllipticalGaussianIntegralWedge(
        double spd,
        double smd,
        double k,
        double x0,
        double y0,
        double r,
#ifdef DEBUG
        int x_sign,
        int y_sign,
#endif
        bool calculate_first_deriv,
        bool calculate_second_deriv,
        double background
    ) :
        EllipticalGaussianIntegralByOrder(calculate_first_deriv,
                                          calculate_second_deriv),
        __wedge_coef(5),
        __integral(1, x0 / r, y0 / r),
        __multipliers(5)
    {
        double x02 = x0 * x0,
               y02 = y0 * y0,
               r2 = r * r,
               x02_p_y02 = x02 + y02;
        if(r2 <= x02_p_y02) {
            __psf_area = 0;
            set_parameters(spd, smd, k, x0, y0, 0, 0, background);
        } else {
            double dx_dy_common = std::sqrt((r2 - x02_p_y02) / x02_p_y02);
            __center_psf_r2 = r2 * std::exp(-0.5 * (spd * x02 + smd * y02)
                                            -
                                            k * x0 * y0);
            double area = __integral(0, 0);
            __psf_area = area * __center_psf_r2;
            set_parameters(spd,
                           smd,
                           k,
                           x0,
                           y0,
                           y0 * dx_dy_common,
                           x0 * dx_dy_common,
                           area * r2 * background);
            if(__psf_area == 0) return;
            for(size_t i = 0; i < 5; ++i) {
                __wedge_coef[i] = __coef[i] * (i < 3 ? r2 : r);
                __multipliers[i].reserve(__initial_storage);
                __multipliers[i].push_back(1);
            }
#ifdef DEBUG
            __dx=y0 * dx_dy_common;
            __dy=x0 * dx_dy_common;
            __radius = r;
            __x_sign = x_sign;
            __y_sign = y_sign;
#endif 
        }
    }

#ifdef DEBUG
    void EllipticalGaussianIntegralWedge::describe(std::ostream &os) const
    {
        os << "wedge (" << std::abs(__x0 - __dx) * __x_sign << ", "
            << std::abs(__y0 - __dy) * __y_sign
            << ", " << __radius << "), orders=";
        for(size_t i = 0; i < __orders.size(); ++i)
            os << __orders[i] << ", ";
        os << "value=" << __value[NO_DERIV] << " +- " << max_error()
            << std::endl;
    }
#endif

} //End PSF namespace.
