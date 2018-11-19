/**\file
 *
 * \brief Defines a basic PSF based on local second order expansion.
 *
 * \ingroup PSF
 */

#ifndef LOCAL_POLYNOMIAL_PSF_H
#define LOCAL_POLYNOMIAL_PSF_H

#include "../Core/SharedLibraryExportMacros.h"
#include "PSF.h"

namespace PSF {

    ///A base class for PSF models which are locally approximated by a polynomial
    ///of up to some degree defined at construction.
    class LIB_PUBLIC LocalPolynomial : public PSF {
    private:
        unsigned
            ///See min_poly_degree argument to constructor.
            __min_poly_degree,

            ///See max_poly_degree argument to constructor.
            __max_poly_degree;

    protected:
        ///\brief Calculates the integral of the PSF over a rectangle.
        ///
        ///Using the local polynomial approximation around the center of the
        ///rectangle.
        virtual double integrate_rectangle(
            ///The x coordinate of the center of the rectangle.
            double center_x,

            ///The y coordinate of the center of the rectangle.
            double center_y,

            ///The full size of the rectangle along x.
            double dx,

            ///The full size of the rectangle along y.
            double dy
        ) const;

        ///\brief Integrates the PSF a wedge of a circle.
        //
        ///The wedge is defined by the following boundaries:
        /// * the line x=x
        /// * the line y=y
        /// * the circle centered at (0, 0) with a radius=radius
        ///If x is 0 the the left vs right wedge is chosen according to left
        ///Same for y0 and bottom.
        virtual double integrate_wedge(
            ///The left/right boundary of the wedge if the given value is
            ///positive/negative.
            double x,

            ///The bottom/top boundary of the wedge if the given value is
            ///positive/negative.
            double y,
            
            ///The radius of the circle giving the rounded boundary of the
            ///wedge.
            double radius, 
            
            ///If x is exactly 0, this values is used to decide to which side of
            ///x the wedge extends.
            bool left=false,

            ///If y is exactly 0, this values is used to decide to which side of
            ///y the wedge extends.
            bool bottom=false
        ) const;

    public:
        ///Create a polynomial approximated PSF of up to the given degree.
        LocalPolynomial(
            ///The maximum order of polynomial coefficients to consider.
            unsigned max_poly_degree,

            ///The minimum order of polynomial coefficients to consider.
            unsigned min_poly_degree=0
        ):
            __min_poly_degree(min_poly_degree),
            __max_poly_degree(max_poly_degree) {}

        ///\brief The coefficient of the term in front of the
        /// \f$x^{x_power}\times y^{y_power}\f$ term 
        ///in the local polynomial approximation of the PSF valid arount x, y
        virtual double poly_coef(
            ///The value of x to evaluate the polynomial term at.
            double x,

            ///The value of y to evaluate the polynomial term at.
            double y,
            
            ///The powerlaw index of x in the term.
            unsigned x_power, 

            ///The powerlaw index of y in the term.
            unsigned y_power
        ) const =0;

        ///Changes the maximum/minimum degree of polynomial terms.
        virtual void set_degree_range(
            ///See same name argument to constructor.
            unsigned max_poly_degree,

            ///See same name argument to constructor.
            unsigned min_poly_degree = 0
        )
        {
            __min_poly_degree=min_poly_degree;
            __max_poly_degree=max_poly_degree;
        }
    };

} //End PSF namespace.

#endif
