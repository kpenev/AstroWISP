/**\file
 *
 * \brief Class approximating integrals of ellptical Gaussians by incremental
 * expansion order.
 *
 * \ingroup PSF
 */

#ifndef __ELLIPTICAL_GAUSSIAN_INTEGRAL_BY_ORDER_H
#define __ELLIPTICAL_GAUSSIAN_INTEGRAL_BY_ORDER_H

#include "../Core/SharedLibraryExportMacros.h"
#include "WedgeIntegral.h"
#include <valarray>

namespace PSF {

    ///\brief Tags for derivatives of integrals of EllipticalGaussianPSF with
    ///respect to the shape parameters.
    enum LIB_PUBLIC SDKDerivative {
        NO_DERIV,	///< The value of the integral
        S_DERIV, 	///< First derivative with respect to S
        D_DERIV,	///< First derivative with respect to D
        K_DERIV,	///< First derivative with respect to K
        SS_DERIV,	///< Second derivative with respect to S
        SD_DERIV,	///< Derivative wrt to D(S) of the S(D) derivative.
        SK_DERIV,	///< Derivative wrt to K(S) of the S(K) derivative.
        DD_DERIV,	///< Second derivative with respect to D
        DK_DERIV,	///< Derivative wrt to K(D) of the D(K) derivative.
        KK_DERIV	///< Second derivative with respect to K
    };

    ///\brief Elliptical gaussian integrals refined by incrementing the
    ///expansion order.
    ///
    ///The area over which integrals are being calculated is assumed to be
    ///small enough for the error bounds in the
    ///<a href="SubPixPhot.pdf">description</a> to apply.
    ///
    /// \ingroup PSF
    class LIB_PUBLIC EllipticalGaussianIntegralByOrder {
    private:
        ///The tightest upper limit to the error in the integral estimate.
        double __error;

        std::valarray<double>
            ///The current S values from the description
            __sum,

            ///The current upper limits on the values of the delta quantities
            __delta;

    protected:
        ///\brief The %PSF at the center of the rectangle times the area of
        ///the rectangle.
        double __psf_area,

               ///\brief The x offset of the point around which the PSF is
               //expanded, relative to the center of the %PSF.
               __x0,

               ///\brief The y offset of the point around which the PSF is
               //expanded, relative to the center of the %PSF.
               __y0
#ifdef DEBUG
                   ,
               ///\brief The half-span of a tight rectangle enclosing the
               ///integral\ area in the x direction.
               __dx,

               ///\brief The half-span of a tight rectangle enclosing the
               ///integral area in the y direction.
               __dy
#endif
                   ;

        ///\brief The values of the integral and its up to second S, D, K
        ///derivatives.
        std::valarray<double> __value;

        ///\brief The orders up to which the terms of equations (59)-(63)
        ///are estiamted
        std::valarray<unsigned> __orders;

        bool __first_deriv, ///< Calculate first order S,D,K derivatives?
             __second_deriv;///< Calculate 2nd order S,D,K derivatives?

        std::valarray<double>
            ///\brief \f$ C_{20}/2 \f$, \f$ C_{11}/2 \f$, \f$ C_{02}/2 \f$,
            /// \f$ C_{10}/2 \f$, \f$ C_{01}/2 \f$
            __coef,

            ///The argument of the exponent of each term of eq 45
            __factors,

            ///The last term of each sum
            __last_term;

        ///\brief Adds what needs to be added to the integral and derivatives
        ///if to_add is being added to the integral.
        ///
        ///Regardless of the region being integrated, updates to the
        ///derivatives depend only on the update to the value of the
        ///integral.
        void fill_new_terms(
                ///The value being added to the estimate of the integral.
                const double to_add,

                ///The (2,0) index of the new term being added.
                unsigned i,

                ///The (1,1) index of the new term being added.
                unsigned j,

                ///The (0,2) index of the new term being added.
                unsigned k,

                ///The (1,0) index of the new term being added.
                unsigned l,

                ///The (0,1) index of the new term being added.
                unsigned m,

                ///The updates to the values of the integral and derivatives,
                ///indexed by SDKDerivative.
                std::valarray<double> &new_terms);

        ///\brief An upper limit to the error in the current integral
        ///estimate.
        ///
        ///Assumes that all members are consistently initialized, sets the
        ///value of __error and returns the tightest error bound.
        double calculate_error();

        ///\brief Use one order higher in the given index estimates of the
        ///integral and any derivatives.
        virtual void update_estimates(unsigned index)=0;

        ///Returns the index of the best term to refine.
        unsigned refine_index();

        ///\brief Increment the order of the index-th term by 1 updating all
        ///members.
        void next_order(unsigned index);
    public:
        ///\brief Default constructor, set_parameters must be called before
        ///first use.
        EllipticalGaussianIntegralByOrder(bool calculate_first_deriv=false,
                bool calculate_second_deriv=false) :
            __first_deriv(calculate_first_deriv),
            __second_deriv(calculate_second_deriv) {}

        ///\brief Sets the parameters of the integral to calculate.
        ///
        ///The protected member __psf_area must already conain the value of
        ///the PSF at the point arounf which the expansion is performed
        ///times the area being integrated.
        void set_parameters(double spd, double smd, double k,
                double x0, double y0, double dx, double dy, double bg_area);

        ///\brief Increments one of I,J,K,L,M
        ///(see <a href="SubPixPhot.pdf">description</a>) by 1.
        ///
        ///Updates the approximations (of the integral and its derivatives)
        ///and errors, returning the new best estimate value of the integral.
        double refine();

        ///\brief Returns the best estimate for the value of the integral
        ///or a derivative.
        double value(SDKDerivative deriv=NO_DERIV) const;

        ///\brief Returns a strict upper limit to the error in the integral
        ///(or a derivative).
        double max_error() const {return __error;}

        ///\brief Returns the orders to which the various sums have been
        ///approximated.
        const std::valarray<unsigned> &orders() const {return __orders;}

#ifdef DEBUG
        ///Output a string describing the integral to a stream.
        virtual void describe(std::ostream &os) const =0;
#endif

        virtual ~EllipticalGaussianIntegralByOrder() {}
    }; //End EllipticalGaussianIntegralByOrder class.

} //End PSF namespace.

#endif
