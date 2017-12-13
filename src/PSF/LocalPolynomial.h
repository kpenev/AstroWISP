/**\file
 *
 * \brief Defines a basic PSF based on local second order expansion.
 *
 * \ingroup PSF
 */

#ifndef LOCAL_POLYNOMIAL_PSF_H
#define LOCAL_POLYNOMIAL_PSF_H

#include "PSF.h"

namespace PSF {

    ///A base class for PSF models which are locally approximated by a polynomial
    ///of up to some degree defined at construction.
    class LocalPolynomial : public PSF {
    private:
        ///The minimum order of polynomial coefficients to consider.
        unsigned __min_poly_degree,

                 ///The maximum order of polynomial coefficients to consider.
                 __max_poly_degree;

    protected:
        ///\brief Calculates the integral of the PSF over a rectangle.
        ///
        ///Using the local polynomial approximation around the center of the
        ///rectangle.
        virtual double integrate_rectangle(double center_x, double center_y,
                double dx, double dy) const;

        ///\brief Integrates the PSF a wedge of a circle.
        //
        ///The wedge is defined by the following boundaries:
        /// * the line x=x
        /// * the line y=y
        /// * the circle centered at (0, 0) with a radius=radius
        ///If x is 0 the the left vs right wedge is chosen according to left
        ///Same for y0 and bottom.
        virtual double integrate_wedge(double x, double y, double radius, 
                bool left=false, bool bottom=false) const;

    public:
        ///Create a polynomial approximated PSF of up to the given degree.
        LocalPolynomial(unsigned max_poly_degree, unsigned min_poly_degree=0):
            __min_poly_degree(min_poly_degree),
            __max_poly_degree(max_poly_degree) {}

        ///The coefficient of the term in front of the x^x_power*y^y_power term 
        ///in the local polynomial approximation of the PSF valid arount x, y
        virtual double poly_coef(double x, double y, unsigned x_power, 
                unsigned y_power) const =0;

        ///Changes the maximum degree of polynomial terms to the given value.
        virtual void set_degree_range(unsigned max_poly_degree,
                                      unsigned min_poly_degree = 0) 
        {
            __min_poly_degree=min_poly_degree;
            __max_poly_degree=max_poly_degree;
        }
    };

} //End PSF namespace.

#endif
