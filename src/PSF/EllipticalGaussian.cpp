/**\file
 *
 * \brief Defining the more complicated members of the elliptical Gaussian
 * point spread function clasess.
 * (%PSF).
 *
 * \ingroup PSF
 */

#include "EllipticalGaussian.h"

namespace PSF {

    double EllipticalGaussian::__default_abs_precision=0,
           EllipticalGaussian::__default_rel_precision=1e-5,
           EllipticalGaussian::__default_max_exp_coef=1;

    Split::Split(double spd, double smd, double k, double x0, double y0,
            double dx, double dy, double max_exp_coef)
    {
        __x_split=std::max(
            1,
            static_cast<int>(
                std::ceil(std::sqrt(spd / 2.0) * dx / max_exp_coef)
            )
        );
        __y_split=std::max(
            1,
            static_cast<int>(
                std::ceil(std::sqrt(smd / 2.0) * dy / max_exp_coef)
            )
        );
        __x_split=std::max(
            __x_split,
            static_cast<int>(
                std::ceil(std::abs(spd * x0 + k * y0) * dx / max_exp_coef)
            )
        );
        __y_split=std::max(
            __y_split,
            static_cast<int>(
                std::ceil(std::abs(smd * y0 + k * x0) * dy / max_exp_coef)
            )
        );
        if(__x_split <= 1 && __y_split <= 1) {
            __sub_x0 = x0;
            __sub_y0 = y0;
            __sub_dx = dx;
            __sub_dy = dy;
            __split = false;
        } else {
            __sub_x0 = x0 - dx / 2.0 * (__x_split - 1) / __x_split;
            __sub_y0 = y0 - dy / 2.0 * (__y_split - 1) / __y_split;
            __sub_dx = dx / __x_split;
            __sub_dy = dy / __y_split;
            __split=true;
        }
    }

    ///Calculates n!
    double factorial(unsigned n)
    {
        double result=1.0;
        for(unsigned i=1; i<=n; i++) result*=i;
        return result;
    }

    void EllipticalGaussian::update_values() const
    {
        __integral_values=__integral_error=0.0;
        for(size_t i=0; i<__pieces.size(); i++) {
            __integral_error+=__pieces[i]->max_error();
            __integral_values[NO_DERIV]+=__pieces[i]->value();
            for(int d=(__first_deriv ? S_DERIV : SS_DERIV);
                    d<=(__second_deriv ? KK_DERIV : K_DERIV); d++)
                __integral_values[d]+=
                    __pieces[i]->value(static_cast<SDKDerivative>(d));
        }
    }

    void EllipticalGaussian::refine() const
    {
        size_t refine_index = 0;
        double max_error=-1.0;
        for( size_t i = 0; i < __pieces.size(); i++ ) {
            if(__pieces[i]->max_error()>max_error) {
                max_error=__pieces[i]->max_error();
                refine_index=i;
            }
        }
        __pieces[refine_index]->refine();
        update_values();
    }

    double EllipticalGaussian::poly_coef(double x, double y,
            unsigned x_power, unsigned y_power) const
    {
        double f0 = operator()(x, y);
        if(x_power == 0 && y_power == 0) return f0;
        f0 -= __bg;
        double c10 = __spd * x + __k * y,
               c01 = __smd * y + __k * x,
               result = 0.0,
               c20_i = 1.0,
               ifact = 1;
        for(unsigned i = 0; i <= x_power / 2; ++i) {
            double c02_k = 1.0,
                   kfact = 1;
            for(unsigned k = 0; k <= y_power/2; ++k) {
                double c11_j = 1.0,
                       c10_l = std::pow(c10,
                                        static_cast<int>(x_power - 2 * i)),
                       c01_m = std::pow(c01,
                                        static_cast<int>(y_power - 2 * k)),
                       jfact = 1.0;
                for(
                    unsigned j = 0;
                    j <= std::min(x_power - 2 * i, y_power - 2 * k);
                    j++
                ) {
                    unsigned l = x_power - 2 * i - j,
                             m = y_power - 2 * k - j;
                    result += (
                        ((i + j + k + l + m) % 2 ? -1 : 1)
                        *
                        c20_i * c11_j * c02_k * c10_l * c01_m
                        /
                        (ifact * jfact * kfact * factorial(l) * factorial(m))
                    );
                    c11_j *= __k;
                    c10_l /= c10;
                    c01_m /= c01;
                    jfact *= j + 1;
                }
                c02_k *= __smd / 2.0;
                kfact *= k + 1;
            }
            c20_i *= __spd / 2.0;
            ifact *= i + 1;
        }
        return f0 * result;
    }

    double EllipticalGaussian::integrate_rectangle(double center_x,
            double center_y, double dx, double dy) const
    {
        Split split(__spd,
                    __smd,
                    __k,
                    center_x,
                    center_y,
                    dx,
                    dy,
                    __max_exp_coef);
        if(split) __pieces.reserve(__pieces.size() + 4 * split.num_pieces());
        for(int yi = 0; yi < split.y_split(); ++yi)
            for(int xi = 0; xi < split.x_split(); ++xi)
                __pieces.push_back(
                        new EllipticalGaussianIntegralRectangle(
                            __spd,
                            __smd,
                            __k,
                            split.sub_x0() + xi * split.sub_dx(),
                            split.sub_y0() + yi * split.sub_dy(),
                            split.sub_dx() / 2,
                            split.sub_dy() / 2,
                            __first_deriv,
                            __second_deriv,
                            __bg
                        )
                );
        return 0;
    }

    double EllipticalGaussian::integrate_wedge(double x_min, double y_min,
            double radius, bool left, bool bottom) const
    {
        int x_sign = (x_min < 0 ? -1 : 1),
            y_sign = (y_min < 0 ? -1 : 1);
        if(x_min == 0 && left) x_sign = -1;
        if(y_min == 0 && bottom) y_sign = -1;

        double r2 = std::pow(radius, 2),
               x_max = x_sign * std::sqrt(r2 - y_min * y_min),
               y_max = y_sign * std::sqrt(r2 - x_min * x_min);
        if(
            std::abs(x_max) <= std::abs(x_min)
            ||
            std::abs(y_max) <= std::abs(y_min)
        )
            return 0;
        double x0 = (x_max + x_min) / 2.0,
               y0 = (y_max + y_min) / 2.0,
               dx = std::abs(x_max - x_min),
               dy = std::abs(y_max - y_min);
        Split split(__spd, __smd, __k, x0, y0, dx, dy, __max_exp_coef);
        if(split) {
            __pieces.reserve(__pieces.size() + 4 * split.num_pieces());
            for(int yi = 0; yi < split.y_split(); ++yi)
                for(int xi = 0; xi < split.x_split(); ++xi)
                    PSF::integrate(split.sub_x0() + xi * split.sub_dx(),
                                   split.sub_y0() + yi * split.sub_dy(),
                                   split.sub_dx(),
                                   split.sub_dy(),
                                   radius);
        } else
            __pieces.push_back(
                new EllipticalGaussianIntegralWedge(__spd,
                                                    __smd,
                                                    __k * (x_sign * y_sign),
                                                    std::abs(x0),
                                                    std::abs(y0),
                                                    radius,
#ifdef DEBUG
                                                    x_sign,
                                                    y_sign,
#endif
                                                    __first_deriv,
                                                    __second_deriv,
                                                    __bg)
            );
        return 0;
    }

    double EllipticalGaussian::find_max(double x1, double y1, double x2,
                    double y2) const
    {
        if(x1 * x2 <= 0 && y1 * y2 <= 0) return 1.0;
        double mk_spd = -__k / __spd,
               mk_smd = -__k / __smd,
               max_x1 = mk_spd * y1,
               max_x2 = mk_spd * y2,
               max_y1 = mk_smd * x1,
               max_y2 = mk_smd * x2,
               max_arg;

        if(max_x1 < x1) max_arg = exp_argument(x1, y1);
        else if(max_x1 > x2) max_arg = exp_argument(x2, y1);
        else max_arg = exp_argument(max_x1, y1);

        if(max_y2 > y2) max_arg = std::max(max_arg, exp_argument(x2, y2));
        else if(max_y2 > y1) max_arg = std::max(max_arg,
                                                exp_argument(x2, max_y2));

        if(max_x2 < x1) max_arg = std::max(max_arg, exp_argument(x1, y2));
        else if(max_x2 < x2) max_arg = std::max(max_arg,
                                                exp_argument(max_x2, y2));

        if(max_y1 > y1 && max_y1 < y2)
            max_arg = std::max(max_arg, exp_argument(x1, max_y1));

        return std::exp(max_arg);
    }

    double EllipticalGaussian::find_min(double x1, double y1, double x2,
                    double y2) const
    {
        return std::exp(std::min(
                    std::min(exp_argument(x1, y1), exp_argument(x1, y2)),
                    std::min(exp_argument(x2, y1), exp_argument(x2, y2))));
    }

    double EllipticalGaussian::integrate(double center_x,
                                            double center_y,
                                            double dx,
                                            double dy,
                                            double circle_radius,
                                            bool calculate_first_deriv,
                                            bool calculate_second_deriv
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                            ,
                                            bool reset_piece_id,
                                            bool skip_piece
#endif
#endif
    ) const
    {
        __first_deriv=calculate_first_deriv;
        __second_deriv=calculate_second_deriv;
        PSF::integrate(center_x, center_y, dx, dy, circle_radius);
        update_values();
        while(__integral_error>std::max(__abs_precision,
                    __integral_values[NO_DERIV]*__rel_precision/
                    (1.0+__rel_precision)))
            refine();
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
        static unsigned piece_id=0;
        if(reset_piece_id) piece_id=0;
#endif
#endif
        for(Core::vector_size_type i=0; i<__pieces.size(); ++i) {
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            if(!skip_piece) {
                std::cerr.precision(16);
                std::cerr.setf(std::ios_base::scientific);
                std::cerr << std::endl << "Piece " << piece_id++ << ": ";
                __pieces[i]->describe(std::cerr);
            }
#endif
#endif
            delete __pieces[i];
        }
        __pieces.clear();
        return __integral_values[NO_DERIV];
    }

}
