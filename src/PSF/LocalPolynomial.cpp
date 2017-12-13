#include "LocalPolynomial.h"

namespace PSF {

    double LocalPolynomial::integrate_rectangle(double center_x,
                                                   double center_y,
                                                   double dx,
                                                   double dy) const
    {
        double integral=0.0, dxh=dx/2, dyh=dy/2, dxh2=dxh*dxh, dyh2=dyh*dyh;
        for(unsigned tot_pow=2*((__min_poly_degree+1)/2);
                tot_pow<=__max_poly_degree; tot_pow+=2) {
            double dxh_p=dxh,
                   dyh_p=std::pow(dyh, static_cast<int>(tot_pow+1));
            for(unsigned x_pow=0; x_pow<=tot_pow; x_pow+=2) {
                double y_pow=tot_pow-x_pow;
                integral+=(4.0*poly_coef(center_x, center_y, x_pow, y_pow)/
                        ((x_pow+1)*(y_pow+1))*dxh_p*dyh_p);
                dxh_p*=dxh2;
                dyh_p/=dyh2;
            }
        }
        return integral;
    }

    double LocalPolynomial::integrate_wedge(double x,
                                               double y,
                                               double radius,
                                               bool left,
                                               bool bottom) const
    {
        short xsign = (x >= 0 ? 1 : -1),
              ysign = (y >= 0 ? 1 : -1);
        if(x == 0 && left) xsign = -1;
        if(y == 0 && bottom) ysign = -1;
        //By symmetry x and y can be assumed positive and we can take care of
        //the other cases by the above signs
        double x0 = std::abs(x),
               y0 = std::abs(y), 
               r2 = radius * radius,
               r4 = r2 * r2,
               x02 = x0 * x0,
               y02 = y0 * y0,
               x03 = x02 * x0,
               y03 = y02 * y0;
        if(x02 + y02 >= r2) return 0.0;
        double xmax2 = r2 - y02,
               xmax = std::sqrt(xmax2),
               xmax3 = xmax2 * xmax,
               ymax2 = r2 - x02,
               ymax = std::sqrt(ymax2),
               ymax3 = ymax2 * ymax,
               gamma = std::atan(xmax / y0) - std::atan(x0 / ymax),
               poly_x = (x0 + xmax) / 2 * xsign,
               poly_y = (y0 + ymax) / 2 * ysign,
               k00 = poly_coef(poly_x, poly_y, 0, 0),
               k10 = poly_coef(poly_x, poly_y, 1, 0) * xsign,
               k01 = poly_coef(poly_x, poly_y, 0, 1) * ysign,
               k11 = poly_coef(poly_x, poly_y, 1, 1) * xsign * ysign,
               k20 = poly_coef(poly_x, poly_y, 2, 0),
               k02 = poly_coef(poly_x, poly_y, 0, 2);
        poly_x *= xsign;
        poly_y *= ysign;
        //Coefficients corrected for the fact that we are not expanding
        //around the corner of the wedge.
        k00 -= (poly_x * k10
                +
                poly_y * k01
                -
                poly_x * poly_x * k20
                -
                poly_y * poly_y * k02
                -
                poly_x * poly_y * k11);
        k10 -= (2.0 * poly_x * k20 + poly_y * k11);
        k01 -= (2.0 * poly_y * k02 + poly_x * k11);
        double result = (
            k00 * 0.5 * (r2 * gamma - y0 * (xmax - x0) - x0 * (ymax - y0))
            +
            k10 * ((ymax3 - y03) / 3.0 - y0 * (xmax2 - x02) / 2.0)
            +
            k01 * ((xmax3 - x03) / 3.0 - x0 * (ymax2 - y02) / 2.0)
            +
            k11 * (r2 * (r2 - x02 - y02)
                   -
                   xmax2 * y02
                   -
                   ymax2 * x02
                   +
                   2.0 * x02 * y02) / 8.0
            +
            k20 * (r4 * gamma / 8.0
                   +
                   y0 * x03 / 3.0
                   -
                   5.0 * y0 * xmax * r2 / 24.0
                   +
                   y03 * xmax / 12.0
                   +
                   x0 * ymax * r2 / 8.0
                   -
                   x03 * ymax / 4.0)
            +
            k02 * (r4 * gamma / 8.0
                   +
                   x0 * y03 / 3.0
                   -
                   5.0 * x0 * ymax * r2 / 24.0
                   +
                   x03 * ymax / 12.0
                   +
                   y0 * xmax * r2 / 8.0
                   -
                   y03 * xmax / 4.0)
        );
        return result;
    }
} //End PSF namespace.
