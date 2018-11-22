/**\file
 *
 * \brief Declares a class that calculates x^my^n integrals over circle
 * pieces.
 *
 * \ingroup PSF
 */

#ifndef __CIRCLE_PIECE_INTEGRAL_H
#define __CIRCLE_PIECE_INTEGRAL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "IntegralUtil.h"
#include "../Core/Typedefs.h"
#include "../Core/Error.h"
#include <cmath>
#include <vector>
#include <cassert>

namespace PSF {

    /**\brief Implements integral of \f$x^my^n\f$ over a horizontal circlre
     * piece in the first quadrant.
     *
     * The wedge is defined by the following diagram:
     *
     * ![](hcircle_piece_diagram.png)
     *
     * The is solution described [here](@ref PSF_integrals_page).
     *
     * \ingroup PSF
     */
    class LIB_LOCAL CirclePieceIntegral {
    private:
        double
            ///Circle center y coordinate (\f$y_c^2\f$) (see diagram).
            __yc,

            ///Square of the circle radius
            __r2,

            ///\brief Square of the circle center y coordinate (\f$y_c^2\f$)
            ///(see diagram)
            __yc2,

            /// \f$r^2-y_c^2\f$
            __r2_m_yc2,

            /// \f$y_{min}-y_c\f$
            __ymin_m_yc,

            /// \f$y_{max}-y_c\f$
            __ymax_m_yc,

            /// \f$\sqrt{r^2-(y_c-y_{min})^2}\f$
            __root_ymin,

            /// \f$\sqrt{r^2-(y_c-y_{max})^2}\f$
            __root_ymax;

        std::vector<double>
            ///Various powers of \f$y_{min}\f$
            __ymin_pow,

            ///Various powers of \f$y_{max}\f$
            __ymax_pow,

            ///Various powers of \f$x_c\f$
            __xc_pow,

            ///Various powers of \f$x_{min}\f$
            __xmin_pow,

            /// \f$R_{n+1}(y_{min})\f$
            __R_ymin,

            /// \f$R_{n+1}(y_{max})\f$
            __R_ymax;

        std::vector< std::vector<double> >
            ///All previously computed values of the integral.
            __values,

            ///\brief All previously computed \f$Q_{m,n}\f$ values.
            ///
            ///The first (outer) index is m and the second (inner) is n. This
            ///way __q[0] is \f$P^{even}\f$ and __q[1] is \f$P^{odd}\f$.
            __q;

        ///\brief Fills \f$R_n(y)\f$.
        void fill_R(
            ///The index (n) up to which to compute \f$R_n(y)\f$ values.
            Core::vector_size_type max_n,

            ///The array of powers of the y argument.
            std::vector<double> &y_pow,

            ///The vector to fill with the computed values, it must
            ///contain at least the \f$R_0(y)\f$ entry.
            std::vector<double> &Rn
        );

        ///\brief Computes all \f$R_n(y_{max})\f$ values for n<=max_n and
        ///stores them in __R_max.
        void fill_R_max(
            ///One past the maximum value of the index for which to compute
            /// \f$R_n\$.
            Core::vector_size_type max_n
        );

        ///\brief Computes all \f$P^{even}_n\f$  for n<max_n and stores them
        ///in __q[0]
        void fill_p_even(
            ///One past the maximum value of teh index for which to compute
            /// \f$P^{even}_n\f$.
            Core::vector_size_type max_n
        );

        ///\brief Computes all \f$P^{odd}_n\f$  for n<max_n and stores tham
        ///in __q[1]
        void fill_p_odd(
            ///One past the maximum value of teh index for which to compute
            /// \f$P^{odd}_n\f$.
            Core::vector_size_type max_n
        );

        ///\brief Fills in a diagonal of __q with indices that sum up to
        ///(m+n) and the first of which has the same parity as m, up to
        ///(m, n), assuming all earlier diagonals are filled.
        void fill_q_diagonal(
            ///See description.
            Core::vector_size_type m,

            ///See description.
            Core::vector_size_type n
        );

        ///\brief Fills in __q[m][n] and possibly other entries of __q if
        ///they are calculated along the way.
        void calculate_q(Core::vector_size_type m, Core::vector_size_type n);

#ifdef DEBUG
        ///\brief A set of asserts verify that the circle piece parameters
        ///define a proper circle piece.
        void assert_parameter_consistency();
#endif

    public:
        ///Create an integral over the area with the given parameters.
        CirclePieceIntegral(
            ///The straight vertical boundary (\f$x_{min}\f$)
            double xmin,

            ///The lower horizontal boundary (\f$y_{min}\f$)
            double ymin,

            ///The upper horizontal boundary (\f$y_{max}\f$)
            double ymax,

            ///The x coordinate of the circle center (\f$x_c\f$)
            double xc,

            ///The y coordinate of the circle center (\f$y_c\f$)
            double yc,

            ///The radius of the circle
            double r,

            ///The integral will be initiall prepared to compute values
            ///for up to this combined x and y power. Calculating values
            ///for larger powers still works, but may result in moving
            ///data around.
            Core::vector_size_type initial_storage = 3
        );

        ///Returns the integral of x^m y^n over the circle piece.
        double operator()(unsigned m, unsigned n);
    }; //End CirclePieceIntegral class.

} //End PSF namespace.

#endif
