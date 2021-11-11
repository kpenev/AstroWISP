/**\file
 *
 * \brief Declares a class that calculates integrals of x^my^n over circle
 * wedges.
 *
 * \ingroup SubPixPhot
 * \ingroup FitSubpix
 */

#ifndef __WEDGE_INTEGRAL_H
#define __WEDGE_INTEGRAL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "IntegralUtil.h"
#include "../Core/Typedefs.h"
#include <cmath>
#include <vector>
#include <cassert>
#include <list>

#ifdef DEBUG
	#include <iostream>
	#include <iomanip>
#endif

namespace PSF {

    /**\brief A class that implements the integral of x^my^n over a circle
     * wedge.
     *
     * The wedge is defined by r, \f$x_0\f$ and \f$y_0\f$ as in the following
     * diagram:
     *
     * ![](circle_wedge.png)
     *
     * The is solution described [here](@ref PSF_integrals_page).
     *
     * \ingroup PSF
     */
    class LIB_LOCAL WedgeIntegral {
    private:
        ///\brief How many entries to allocate storage for in the __values and
        ///__q attributes at construction.
        static const Core::vector_size_type __initial_storage = 30;

        double
            ///\brief The x distance from the circle center to the center of
            ///the chord.
            __x0,

            ///\brief The y distance from the circle center to the center
            ///of the chord.
            __y0,

            __y02,    	///< The square of __y0.
            __twice_y0,	///< 2*__y0
            __deltax,	///< \f$\Delta x\f$
            __deltay,	///< \f$\Delta y\f$
            __r2,    	///< The square of the circle radius.
            __r2_m_y02,	///< \f$r^2-y_0^2\f$.
            __x_min, 	///< The smallest distance to the circle center in x.
            __y_min, 	///< The smallest distance to the circle center in y.
            __x_max, 	///< The largest distance to the circle center in x.
            __y_max, 	///< The largest distance to the circle center in y.

            ///\brief The n independent part of the term from the recursive
            ///relation for \f$P^odd\f$ which does not depend on previous
            /// \f$P^odd\f$ for even n.
            __p_odd_free_term_mult_even,

            ///\brief The n independent part of the term from the recursive
            ///relation for \f$P^odd\f$ which does not depend on previous
            /// \f$P^odd\f$ for odd n.
            __p_odd_free_term_mult_odd;

        std::vector< std::vector<double> >
            ///All previously computed values of the integral.
            __values,

            ///\brief All previously computed \f$Q_{m,n}\f$ values.
            ///
            ///The first (outer) index is m and the second (inner) is n. This
            ///way __q[0] is \f$P^{even}\f$ and __q[1] is \f$P^{odd}\f$.
            __q;

        std::vector<double>
            ///Various powers of \f$\Delta x\f$.
            __deltax_pow,

            ///Various powers of \f$\Delta y\f$.
            __deltay_pow,

            ///Previously computed powers of __y0.
            __x0_pow;

        ///\brief Computes all \f$P^{even}_n\f$  for n<max_n and stores them
        ///in __q[0].
        void fill_p_even(
            ///The index of the last term of \f$P^{even}$ to calculate.
            Core::vector_size_type max_n
        );

        ///\brief Computes all \f$P^{odd}_n\f$  for n<max_n and stores them
        ///in __q[1]
        void fill_p_odd(
            ///The index of the last term of \f$P^{odd}$ to calculate.
            Core::vector_size_type max_n
        );

        ///\brief Fills in a diagonal of __q with indices that sum up to
        ///(m+n) and the first of which has the same parity as m, up to
        ///(m, n), assuming all earlier diagonals are filled.
        void fill_q_diagonal(
            ///The first index of the final \f$Q_{m,n}\f$ term to fill.
            Core::vector_size_type m,

            ///The second index of the final \f$Q_{m,n}\f$ term to fill.
            Core::vector_size_type n
        );

        ///\brief Fills in __q[m][n] and possibly other entries of __q if they
        ///are calculated along the way.
        void calculate_q(
            ///The first index of \f$Q_{m,n}\f$ of the new term to calculate.
            Core::vector_size_type m,

            ///The second index of \f$Q_{m,n}\f$ of the new term to calculate.
            Core::vector_size_type n
        );

#ifdef DEBUG
        ///Outputs the currently calculated \f$Q_{m,n}\f$ values.
        void output_q(
            ///The stream to output to.
            std::ostream &os
        );
#endif

        ///\brief Increases the size of all vectors holding pre-calculated
        ///values by a power of 2 such that they can hold at least the given
        ///m and n values.
        void expand_storage(
            ///The first index of the term for which storage must be guaranteed.
            Core::vector_size_type min_m_size,

            ///The second index of the term for which storage must be
            ///guaranteed.
            Core::vector_size_type min_n_size
        );
    public:
        ///\brief Create an integral over a wedge with the center of the
        ///chord at (x0, y0) in a circle with radius r.
        WedgeIntegral(
            ///The radius of the circle defining the outer boundary of the
            ///wedge.
            double r,

            ///The x coordinate of the center of the cord defining the wedge.
            double x0,

            ///The y coordinate of the center of the cord defining the wedge.
            double y0
        );

        ///Returns the integral of x^m y^n over the wedge.
        double operator()(
            ///The powerlaw index of x in the term to integrate.
            unsigned m,

            ///The powerlaw index of y in the term to integrate.
            unsigned n
        );
    }; //End WedgeIntegral class.

} //Eend PSF namespace.
#endif
