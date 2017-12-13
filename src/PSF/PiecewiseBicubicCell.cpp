/**\file
 * \brief The implementation details of PiecewiseBicubicCell.
 */

#include "PiecewiseBicubicCell.h"
#include <iostream>

namespace PSF {

    std::valarray<double> PiecewiseBicubicCell::integrate_rectangle(
            double xmin, double xmax, double ymin, double ymax,
            const std::valarray<double> &coef_sets) const
    {
        const std::valarray<double> *coef=&coef_sets;
        if(coef_sets.size()==0) coef=&__coef;
        size_t coef_ind=0, num_coef_sets=coef->size()/16;
#ifdef DEBUG
        assert(coef->size()>=16);
        assert(coef->size()%16==0);
#endif
        std::valarray<double> result(0.0, num_coef_sets);
        double ymin_pow=ymin, ymax_pow=ymax;
        for(unsigned ypow=0; ypow<4; ++ypow) {
            double xmin_pow=xmin, xmax_pow=xmax,
                   yfactor=(ymax_pow-ymin_pow)/(ypow+1);
            for(unsigned xpow=0; xpow<4; ++xpow) {
                for(size_t set_ind=0; set_ind<num_coef_sets; ++set_ind)
                    result[set_ind] += (*coef)[16 * set_ind + coef_ind]
                                       *
                                       (xmax_pow - xmin_pow)
                                       /
                                       (xpow + 1)
                                       *
                                       yfactor;
                ++coef_ind;
                xmin_pow*=xmin;
                xmax_pow*=xmax;
            }
            ymin_pow*=ymin;
            ymax_pow*=ymax;
        }
        return result;
    }

    double PiecewiseBicubicCell::integrate_hcircle_piece(double ymin,
            double ymax, double xbound, double radius, double circle_x, 
            double circle_y) const
    {
        bool switch_sign=false;
        if(xbound<circle_x) {
            xbound=-xbound;
            circle_x=-circle_x;
            switch_sign=true;
        }
        if(ymin<0) ymin=0;
        if(ymax>vertical_size()) ymax=vertical_size();
        if(xbound>horizontal_size()) xbound=horizontal_size();
        CirclePieceIntegral integral(xbound, ymin, ymax, circle_x, circle_y,
                radius);
        double result=0;
        unsigned coef_ind=0;
        for(unsigned ypow=0; ypow<=3; ++ypow) {
            int sign=1;
            for(unsigned xpow=0; xpow<=3; ++xpow) {
                result+=__coef[coef_ind]*integral(xpow, ypow)*sign;
#ifdef DEBUG
                assert(!std::isnan(result));
#endif
                if(switch_sign) sign=-sign;
                ++coef_ind;
            }
        }
        return result;
    }

    double PiecewiseBicubicCell::integrate_vcircle_piece(double xmin,
            double xmax, double ybound, double radius, double circle_x,
            double circle_y) const
    {
        bool switch_sign=false;
        if(ybound<circle_y) {
            ybound=-ybound;
            circle_y=-circle_y;
            switch_sign=true;
        }
        if(xmin<0) xmin=0;
        if(xmax>horizontal_size()) xmax=horizontal_size();
        if(ybound>vertical_size()) ybound=vertical_size();
        CirclePieceIntegral integral(ybound, xmin, xmax, circle_y, circle_x,
                radius);
        double result=0;
        unsigned coef_ind=0;
        int sign=1;
        for(unsigned ypow=0; ypow<=3; ++ypow) {
            for(unsigned xpow=0; xpow<=3; ++xpow) {
                result+=__coef[coef_ind]*integral(ypow, xpow)*sign;
                ++coef_ind;
            }
            if(switch_sign) sign=-sign;
        }
        return result;
    }

    std::valarray<double> PiecewiseBicubicCell::operator()(
            double x,
            double y,
            const std::valarray<double> &coef_sets
    ) const
    {
        const std::valarray<double> *coef=&coef_sets;
        if(coef_sets.size()==0) coef=&__coef;

        size_t coef_ind=0, num_coef_sets=coef->size() / 16;
#ifdef DEBUG
        assert(coef->size() >= 16);
        assert(coef->size() % 16 == 0);
#endif

        double y_pow=1;
        std::valarray<double> result(0.0, num_coef_sets);
        for(unsigned ypow=0; ypow<4; ++ypow) {
            double x_pow=1;
            for(unsigned xpow=0; xpow<4; ++xpow) {
                for(size_t set_ind=0; set_ind<num_coef_sets; ++set_ind)
                    result[set_ind] += (*coef)[16 * set_ind + coef_ind]
                                       *
                                       x_pow
                                       *
                                       y_pow;
                ++coef_ind;
                x_pow *= x;
            }
            y_pow *= y;
        }
        return result;
    }

} //End PSF namespace.
