/**\file
 *
 * \brief Defines the methods of the
 * EllipticalGaussianIntegralRectangle class.
 *
 * \ingroup PSF
 */

#include "EllipticalGaussianIntegralRectangle.h"

namespace PSF {

    void EllipticalGaussianIntegralRectangle::update_estimates(
        unsigned index
    )
    {
        double i_term = 1.0;
        const double f3_sq = __factors[3] * __factors[3];
        const double f4_sq = __factors[4] * __factors[4];

        std::valarray<double> new_terms( 0.0, KK_DERIV + 1 );

        double new_term_factor = __psf_area * __last_term[index];

        unsigned i_endindex = __orders[0];
        unsigned i_startindex = index ? 0 : i_endindex;
        for( unsigned i = i_startindex; i <= i_endindex; ++i ) {
            double j_term = i_term;
            unsigned j_endindex = __orders[1];
            unsigned j_startindex = index != 1 ? 0 : j_endindex;

            for( unsigned j = j_startindex; j <= j_endindex; ++j ) {
                double k_term = j_term;
                unsigned k_endindex = __orders[2];
                unsigned k_startindex = index != 2 ? 0 : k_endindex;
                unsigned j_parity = j%2;

                for( unsigned k = k_startindex; k <= k_endindex; ++k ) {
                    double l_term = k_term;

                    if( index != 3 && j_parity ) {
                        l_term *= -__factors[3];
                    }

                    unsigned l_endindex = __orders[3];
                    unsigned l_startindex =
                            index != 3
                        ?   j_parity
                        :   l_endindex + ( l_endindex + j )%2
                        ;

                    for(unsigned l = l_startindex; l <= l_endindex; l += 2) {
                        double m_term = l_term / ( 2 * i + j + l + 1);

                        if ( index!=4 && j_parity ) {
                            m_term *= -__factors[4];
                        }

                        unsigned m_endindex = __orders[4];
                        unsigned m_startindex =
                                index != 4
                            ?   j_parity
                            :   m_endindex + ( m_endindex + j )%2
                            ;

                        for(unsigned m=m_startindex; m<=m_endindex; m+=2 ) {
                            fill_new_terms(m_term/(2*k+j+m+1), i, j, k, l, m,
                                    new_terms);
                            if (m!=m_endindex) {
                                m_term *= f4_sq / ( (m + 1 ) * ( m + 2 ) );
                            }
                        }

                        if ( l != l_endindex ) {
                            l_term *= f3_sq / ( ( l + 1 ) * ( l + 2 ) );
                        }
                    }

                    if ( k != k_endindex ) {
                        k_term *= -__factors[2] / ( k + 1 );
                    }
                }

                if ( j != j_endindex ) {
                    j_term *= -__factors[1] / ( j + 1 );
                }
            }

            if ( i != i_endindex ) {
                i_term *= -__factors[0] / ( i + 1 );
            }
        }

        __value[NO_DERIV] += new_term_factor * new_terms[NO_DERIV];

        int d_startindex = __first_deriv ? S_DERIV : SS_DERIV;
        int d_endindex = __second_deriv ? KK_DERIV : K_DERIV;
        for( int d = d_startindex; d <= d_endindex; ++d ) {
            __value[d] += new_term_factor * new_terms[d];
        }
    }

    EllipticalGaussianIntegralRectangle::EllipticalGaussianIntegralRectangle(
        double spd,
        double smd,
        double k,
        double x0,
        double y0,
        double dx,
        double dy,
        bool calculate_first_deriv,
        bool calculate_second_deriv,
        double background
    ) :
        EllipticalGaussianIntegralByOrder(calculate_first_deriv,
                                          calculate_second_deriv)
    {
        double area=4.0*dx*dy;
        __psf_area=area*std::exp(-0.5*(spd*x0*x0 + smd*y0*y0)-k*x0*y0);
        set_parameters(spd, smd, k, x0, y0, dx, dy, background*area);
#ifdef DEBUG
        __dx=dx;
        __dy=dy;
#endif
    }

#ifdef DEBUG
    void EllipticalGaussianIntegralRectangle::describe(
        std::ostream &os
    ) const
    {
        os << "rectangle (" << __x0-__dx << ", " << __y0-__dy << ") - ("
            << __x0+__dx << ", " << __y0+__dy << "), orders=";
        for(size_t i=0; i<__orders.size(); ++i)
            os << __orders[i] << ", ";
        os << "value=" << __value[NO_DERIV] << " +- " << max_error()
            << std::endl;
    }
#endif

} //End PSF namespace.
