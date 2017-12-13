/**\file
 *
 * \brief The definitions of the methods of the
 * EllipticalGaussianIntegralByOrder class.
 *
 * \ingroup PSF
 */

#include "EllipticalGaussianIntegralByOrder.h"
#include "../Core/NaN.h"

namespace PSF {

    void EllipticalGaussianIntegralByOrder::fill_new_terms(
            const double to_add, 
            unsigned i, unsigned j, unsigned k, unsigned l, unsigned m,
            std::valarray<double> &new_terms)
    {
        new_terms[NO_DERIV] += to_add;

        if(__first_deriv || __second_deriv) {
            const double x0_c3 = __x0 / __coef[3];
            const double x0_c4 = __x0 / __coef[4];
            const double y0_c3 = __y0 / __coef[3];
            const double y0_c4 = __y0 / __coef[4];
            const double two_c0 = 2.0 * __coef[0];
            const double two_c2 = 2.0 * __coef[2];

            double dS_mult=(i ? i/two_c0 : 0),
                   dD_mult = dS_mult,
                   dK_mult = 0;

            if(k) {
                double k_per_two_c2=k/two_c2;
                dS_mult+=k_per_two_c2;
                dD_mult-=k_per_two_c2;
            }
            if(l) {
                double x0_c3_times_l=x0_c3*l;
                dS_mult+=x0_c3_times_l;
                dD_mult+=x0_c3_times_l;
                dK_mult+=y0_c3*l;
            }
            if (m) {
                double y0_c4_times_m=y0_c4*m;
                dS_mult+=y0_c4_times_m;
                dD_mult-=y0_c4_times_m;
                dK_mult+=x0_c4*m;
            }
            if (j) dK_mult+=double(j)/__coef[1];

            if(__first_deriv) {
                new_terms[S_DERIV]+=to_add*dS_mult;
                new_terms[D_DERIV]+=to_add*dD_mult;
                new_terms[K_DERIV]+=to_add*dK_mult;
            }

            if(__second_deriv) {
                double two_c0_pow_2=two_c0*two_c0;
                double two_c2_pow_2=two_c2*two_c2;

                double d2SD_extra = (i ? i/two_c0_pow_2 : 0)
                                    + (k ? k/two_c2_pow_2 : 0);

                if(__x0) {
                    if(l) {
                        const double x0_c3_pow_2=x0_c3*x0_c3;
                        d2SD_extra+=x0_c3_pow_2*l;
                        new_terms[SD_DERIV]-=to_add*x0_c3_pow_2*l;
                    }
                    if(m) new_terms[KK_DERIV]-=to_add*std::pow(x0_c4, 2)*m;
                }

                if(__y0) {
                    if(m) {
                        const double y0_c4_pow_2=y0_c4*y0_c4;

                        d2SD_extra+=y0_c4_pow_2*m;

                        new_terms[SD_DERIV]+=to_add*y0_c4_pow_2*m;
                    }
                    if(l) new_terms[KK_DERIV]-=to_add*std::pow(y0_c3, 2)*l;
                }

                new_terms[SS_DERIV]+=to_add*(std::pow(dS_mult, 2)
                                             -
                                             d2SD_extra);
                new_terms[SD_DERIV]+=to_add*(dS_mult*dD_mult
                                             - (i ? i/two_c0_pow_2 : 0)
                                             + (k ? k/two_c2_pow_2 : 0));
                const double to_add_dK=to_add*dK_mult;

                new_terms[SK_DERIV]+=to_add_dK*dS_mult;
                new_terms[DD_DERIV]+=to_add*(std::pow(dD_mult, 2)
                                             -
                                             d2SD_extra);
                new_terms[DK_DERIV] += to_add_dK * dD_mult;
                new_terms[KK_DERIV] += to_add_dK * dK_mult;

                if(j) new_terms[KK_DERIV]-=to_add*j/std::pow(__coef[1], 2);

                if(__x0 && __y0) {
                    double l_c3_sq=(l ? l/std::pow(__coef[3], 2) : 0),
                           m_c4_sq=(m ? m/std::pow(__coef[4], 2) : 0);

                    new_terms[SK_DERIV]-=to_add*__x0*__y0*(l_c3_sq
                                                           +
                                                           m_c4_sq);
                    new_terms[DK_DERIV]-=to_add*__x0*__y0*(l_c3_sq
                                                           -
                                                           m_c4_sq);
                }
            }
        }
    }

    unsigned EllipticalGaussianIntegralByOrder::refine_index()
    {
        unsigned result = 0;
        double max_term=-1;
        for(unsigned i=0; i<5; i++) {
            double term=1.0;
            for(unsigned j=0; j<5; j++) term*=(i==j ? std::abs(__delta[j]) :
                    std::abs(__sum[j]));
            if(term>max_term) {
                max_term=term;
                result=i;
            }
        }
        return result;
    }

    double EllipticalGaussianIntegralByOrder::calculate_error()
    {
        double delta_plus_s=1.0, s=1.0;
        for(unsigned i=0; i<5; i++) {
            delta_plus_s*=(__delta[i]+std::abs(__sum[i]));
            s*=std::abs(__sum[i]);
        }
        __error=(delta_plus_s-s)*__psf_area;

        return __error;
    }

    void EllipticalGaussianIntegralByOrder::next_order(unsigned index)
    {
        ++__orders[index];
        __last_term[index]*=-__factors[index]/__orders[index];
        if(index==1 || index>2) __sum[index]+=std::abs(__last_term[index]);
        __delta[index]*=std::abs(__factors[index])/(__orders[index]+1);
        update_estimates(index);
        calculate_error();
    }

    void EllipticalGaussianIntegralByOrder::set_parameters(
        double spd,
        double smd,
        double k,
        double x0,
        double y0,
        double dx,
        double dy,
        double bg_area
    )
    {
        __value.resize(KK_DERIV+1, 0);
        if(__psf_area) {
            __sum.resize(5, 1);
            __delta.resize(5);
            __x0=x0;
            __y0=y0;
            __orders.resize(5, static_cast<unsigned>(0));
            __coef.resize(5);
            __factors.resize(5);
            __last_term.resize(5, 1.0);

            __value[NO_DERIV]=__psf_area+bg_area;
            __coef[0]=spd/2.0; __coef[1]=k; __coef[2]=smd/2.0;
            __coef[3]=spd*x0+k*y0; __coef[4]=smd*y0+k*x0;

            __factors[0]=__coef[0]*dx*dx;
            __factors[1]=__coef[1]*dx*dy;
            __factors[2]=__coef[2]*dy*dy;
            __factors[3]=__coef[3]*dx;
            __factors[4]=__coef[4]*dy;

            for(unsigned i=0; i<5; i++)
                __delta[i]=std::abs(__factors[i])*std::exp(
                        (i==0 || i==2
                         ? -__factors[i]
                         : std::abs(__factors[i])));
            __error=Core::Inf;
        } else {
            __error=0;
        }
    }

    double EllipticalGaussianIntegralByOrder::refine()
    {
        while(static_cast<int>(__orders[0])<=
                static_cast<int>(std::ceil(__factors[0]/2-1))) next_order(0);
        while(static_cast<int>(__orders[2])<=
                static_cast<int>(std::ceil(__factors[2]/2-1))) next_order(2);
        next_order(refine_index());
        return __value[NO_DERIV];
    }

    double EllipticalGaussianIntegralByOrder::value(
        SDKDerivative deriv
    ) const
    {
        double x2py2 = (__x0 * __x0 + __y0 * __y0) / 2.0,
               x2my2 = (__x0 * __x0 - __y0 * __y0) / 2.0,
               xy = __x0 * __y0;

        switch(deriv) {
            case NO_DERIV : return __value[NO_DERIV]; 
            case S_DERIV : return (__value[S_DERIV]
                                   -
                                   x2py2 * __value[NO_DERIV]); 
            case D_DERIV : return (__value[D_DERIV]
                                   -
                                   x2my2 * __value[NO_DERIV]);
            case K_DERIV : return __value[K_DERIV] - xy * __value[NO_DERIV];
            case SS_DERIV : return (__value[SS_DERIV]
                                    -
                                    2.0 * x2py2 * __value[S_DERIV]
                                    +
                                    x2py2 * x2py2 * __value[NO_DERIV]);
            case SD_DERIV : return (__value[SD_DERIV]
                                    -
                                    x2py2 * __value[D_DERIV]
                                    -
                                    x2my2 * __value[S_DERIV]
                                    +
                                    x2py2 * x2my2 * __value[NO_DERIV]);
            case SK_DERIV : return (__value[SK_DERIV]
                                    -
                                    x2py2 * __value[K_DERIV]
                                    -
                                    xy * __value[S_DERIV]
                                    +
                                    xy * x2py2 * __value[NO_DERIV]);
            case DD_DERIV : return (__value[DD_DERIV]
                                    -
                                    2.0 * x2my2 * __value[D_DERIV]
                                    +
                                    x2my2 * x2my2 * __value[NO_DERIV]);
            case DK_DERIV : return (__value[DK_DERIV]
                                    -
                                    x2my2 * __value[K_DERIV]
                                    -
                                    xy * __value[D_DERIV]
                                    +
                                    xy * x2my2 * __value[NO_DERIV]);
            case KK_DERIV : return (__value[KK_DERIV]
                                    -
                                    2.0 * xy * __value[K_DERIV]
                                    +
                                    xy * xy * __value[NO_DERIV]);
            default : return Core::NaN;
        }
    }

} //End PSF namespace.
