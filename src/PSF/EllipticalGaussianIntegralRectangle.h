/**\file
 * 
 * \brief Declare a class for calculating ellptical Gaussian integrals over
 * rectangles.
 *
 * \ingroup PSF
 */

#ifndef __ELLIPTICAL_GAUSSIAN_INTEGRAL_RECTANGLE_H
#define __ELLIPTICAL_GAUSSIAN_INTEGRAL_RECTANGLE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EllipticalGaussianIntegralByOrder.h"

namespace PSF {

    ///\brief Calculates integrals of an elliptical gaussion over a 
    ///rectangle.
    ///
    ///\ingroup PSF
    class LIB_PUBLIC EllipticalGaussianIntegralRectangle :
        public EllipticalGaussianIntegralByOrder {
    protected:
        ///\brief Use one order higher in the given index estimates of the
        ///integral and any derivatives.
        void update_estimates(unsigned index);
    public:
        ///\brief Create an integral with the given S+D, S-D and K over a
        ///rectangle.
        ///
        ///The rectangle to integrate over is given by x0-dx < x < x0+dx,
        ///y0-dy < y < y0+dy. Optionally the first (if calculate_first_deriv 
        ///is true) and second (if calculate_second_deriv is true) 
        ///derivatives of the integral are also calculated.
        EllipticalGaussianIntegralRectangle(double spd, double smd, double k,
                double x0, double y0, double dx, double dy,
                bool calculate_first_deriv=false,
                bool calculate_second_deriv=false, double background=0);

#ifdef DEBUG
        ///Output a string describing the integral to a stream. 
        void describe(std::ostream &os) const;
#endif
    }; //End EllipticalGaussianIntegralRectangle class.

} //End PSF namespace.

#endif
