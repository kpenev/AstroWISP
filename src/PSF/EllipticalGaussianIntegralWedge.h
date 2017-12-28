/**\file
 *
 * \brief Declare a class for integrating elliptical Guassians over wedge
 * areas.
 *
 * \ingroup PSF
 */

#ifndef __ELLIPTICAL_GAUSSIAN_INTEGRAL_WEDGE_H
#define __ELLIPTICAL_GAUSSIAN_INTEGRAL_WEDGE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EllipticalGaussianIntegralByOrder.h"

namespace PSF {

    ///\brief Calculates integrals of an elliptical gaussion over a wedge of
    ///a circle.
    ///
    ///\ingroup PSF
    class LIB_PUBLIC EllipticalGaussianIntegralWedge :
        public EllipticalGaussianIntegralByOrder {
    private:
        ///The value of the PSF at the midpoint of the wedge chord.
        double __center_psf_r2;

        std::valarray<double>
            ///\brief \f$ C_{20}r^2/2 \f$, \f$ C_{11}r^2/2 \f$,
            /// \f$ C_{02}r^2/2 \f$, \f$ C_{10}r/2 \f$, \f$ C_{01}r/2 \f$
            __wedge_coef;

#ifdef DEBUG
        double __radius;///< Tha radius of the circle

        int 
            /// +1 for right facing wedges, -1 for left facing ones.
            __x_sign,

            ///+1 for top facing wedges, -1 for bottom facing ones.
            __y_sign;
#endif

        ///\brief Initial storage to reserve for storing quantities required
        ///multiple times.
        static const Core::vector_size_type __initial_storage=100;

        ///\brief The integral being refined until the desired precision is 
        ///achieved.
        WedgeIntegral __integral;

        ///The i,j-th element is (-__coef[i])^j/j!
        std::vector< std::vector<double> > __multipliers;
    protected:
        ///\brief Use one order higher in the given index estimates of the
        ///integral and any derivatives.
        void update_estimates(unsigned index);
    public:
        ///Construct a wedge integral for a PSF with the given parameters.
        EllipticalGaussianIntegralWedge(
                double spd,///< S+D
                double smd,///< S-D
                double k,///< K

                ///The x coordinate of the midpoint of the wedge chord
                ///relative to the circle center.
                double x0,

                ///The y coordinate of the midpoint of the wedge chord
                ///relative to the circle center.
                double y0,

                double r,///< Radius of the circle

#ifdef DEBUG
                ///+1 for right facing wedges, -1 for left facing ones.
                int x_sign,

                ///+1 for top facing wedges, -1 for bottom facing ones.
                int y_sign,
#endif

                ///Should first order derivatives be calculated?
                bool calculate_first_deriv=false,

                ///Should second order derivatives be calculated?
                bool calculate_second_deriv=false,

                ///The background under the PSF.
                double background=0);

#ifdef DEBUG
        ///Output a string describing the integral to a stream. 
        void describe(std::ostream &os) const;
#endif
    };

} //End PSF namespace.

#endif
