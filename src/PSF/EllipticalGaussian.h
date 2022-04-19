/**\file
 *
 * \brief Class defining ellptical Gaussian point spread functions (%PSF).
 *
 * The elliptical gaussian %PSF is such that the intensity of light from an
 * isolated source falling on the detector at a location (x,y) relative to the
 * center of light for the source is given by:
 *
 * \f$ I(x,y) = B + A\exp\left\{ -\frac{1}{2} \left[S(x^2 + y^2) +
 * D(x^2-y^2) + 2Kxy\right]\right\}\f$
 *
 * where S, D and K are shape parameters, B is the background under the
 * source and A is the amplitude (how bright) of the source. For more details
 * see: <a href="http://arxiv.org/abs/0906.3486"> Pal, A.\ 2009, Ph.D.~Thesis
 * </a>
 *
 * Details of how the integration is performed with up to a prescribed
 * precision are given <a href="SubPixPhot.pdf">here</a>.
 *
 * \ingroup PSF
 */

#ifndef __ELLIPTICAL_GAUSSIAN_PSF_H
#define __ELLIPTICAL_GAUSSIAN_PSF_H

#include "../Core/SharedLibraryExportMacros.h"
#include "EllipticalGaussianIntegralRectangle.h"
#include "EllipticalGaussianIntegralWedge.h"
#include "PSF.h"
#include "Split.h"
#include "../Core/NaN.h"
#include <vector>
#include <iostream>

namespace PSF {

    ///\brief An ellptical Gaussian PSF with peak flux above the  background
    ///normalized to 1.
    ///
    /// \ingroup PSF
    class LIB_PUBLIC EllipticalGaussian : public PSF {
    private:
      ///Changed these to be an independent declaration each
        static double
            ///Default for __abs_precision
            __default_abs_precision;
        static double
            ///Default for __rel_precision
            __default_rel_precision;
        static double
            ///Default for __max_exp_coef
            __default_max_exp_coef;

        double __spd,	///< S+D
               __smd,	///< S-D
               __k, 	///< K
               __bg,	///< background (divided by the actual amplitude).

               ///\brief The maximum value individual exponent arguments are
               ///allowed to take.
               ///
               ///Higher values result in faster calculation of integrals,
               ///but worse precision. Should always be of order unity.
               __max_exp_coef;

        mutable double
            __abs_precision,///< The absolute precision required
            __rel_precision, ///< The relative precision required

            ///\brief The tightest upper limit to the error in the integral
            ///estimation.
            __integral_error;

        ///The integrals over individiual pieces of the integration area.
        mutable std::vector<EllipticalGaussianIntegralByOrder*> __pieces;

        ///\brief The current best estimate for the integral and its S,D,K
        ///derivatives.
        mutable std::valarray<double> __integral_values;

        ///Whether to calculate first order derivatives of integrals.
        mutable bool __first_deriv,

             ///Whether to calculate second order derivatives of integrals.
             __second_deriv;

        ///Returns the argument inside the exponent of the PSF.
        inline double exp_argument(double x, double y) const
        {return -0.5*(__spd*x*x+__smd*y*y+2.0*__k*x*y);}

        ///\brief Updates __integral_values and __integral_error according to
        ///the latest state of __pieces
        void update_values() const;

        ///\brief Increments one of I,J,K,L,M
        ///(see <a href="SubPixPhot.pdf">description</a>) by 1 for a single
        ///integration piece.
        ///
        ///Updates the approximations (of the integral and its derivatives)
        ///and errors.
        void refine() const;

        ///\brief Returns the maximum value the PSF minus background takes
        ///over a rectangle.
        ///
        ///The rectangle is between x1 and x2 in the horizontal direction and
        ///between y1 and y2 in the vertical.
        double find_max(
                double x1, double y1, double x2, double y2) const;

        ///\brief Returns the minimum value the PSF minus background takes
        ///over a rectangle.
        ///
        ///The rectangle is between x1 and x2 in the horizontal direction and
        ///between y1 and y2 in the vertical.
        double find_min(double x1, double y1, double x2, double y2) const;
    protected:
        ///\brief Calculates the integral of the PSF over a rectangle.
        ///
        ///Using the local polynomial approximation around the center of the
        ///rectangle. The value of the integral is refined by icreasing the
        ///expansion order.
        virtual double integrate_rectangle(double center_x, double center_y,
                double dx, double dy) const;

        ///Integrates the PSF over a circle wedge.
        ///
        ///The wedge is defined by the following boundaries:
        /// * the line x=x
        /// * the line y=y
        /// * the circle centered at (0, 0) with a radius=radius
        ///If x is 0 the the left vs right wedge is chosen according to left
        ///Same for y0 and bottom.
        ///
        ///The value of the integral is refined by increasing the expansion
        ///order, rather than by subdivisions.
        virtual double integrate_wedge(double x_min,
                                       double y_min,
                                       double radius,
                                       bool left = false,
                                       bool bottom = false) const;
    public:
        ///\brief Construct an elliptical Gaussian PSF with the given
        ///parameters.
        ///
        ///Integrals are calculated to the less stringent of the absolute or
        ///relative precisions given.
        EllipticalGaussian(
                ///The coefficient in front of \f$(x^2+y^2)/2\f$ in the
                ///exponent
                double s,

                ///The coefficient in front of \f$(x^2-y^2)/2\f$ in the
                ///exponent
                double d,

                ///The coefficient in front of \f$xy\f$ in the exponent.
                double k,

                ///Background to add under the PSF.
                double background=0,

                ///The relative precision to impose on integrals.
                double relative_precision=__default_rel_precision,

                ///The absolute precision to impose on integrals.
                double absolute_precision=__default_abs_precision,

                ///The maximum value to allow for any term appearing in the
                ///exponent when evaluating integrals. Larger values result in
                ///faster code, but very large values lead to severe numerical
                ///roundoff errors. Values of order unity work well.
                double max_exp_coef=__default_max_exp_coef):
            __spd(s+d), __smd(s-d), __k(k), __bg(background),
            __max_exp_coef(max_exp_coef), __abs_precision(absolute_precision),
            __rel_precision(relative_precision),
            __integral_values(KK_DERIV+1) {}

        ///\brief Change the precision requirements from what was set at
        ///construction.
        void set_precision(double relative, double absolute) const
        {__abs_precision=absolute; __rel_precision=relative;}

        ///Change the default precision requirements.
        static void set_default_precision(double relative, double absolute)
        {__default_abs_precision=absolute; __default_rel_precision=relative;}

        ///\brief Change the default max exponent arguments allowed during
        ///calculations.
        static void set_default_max_exp_coef(double value)
        {__default_max_exp_coef=value;}

        ///Evaluates the PSF at the given position
        inline double operator()(double x, double y) const
        {double result; evaluate(x, y, &result); return result;}

        ///Calculates the PSF and its derivatives at the given position
        template<class ITERABLE>
        void evaluate(
                ///The x coordinate where to evaluate the PSF (and
                ///derivatives).
                double x,

                ///The y coordinate where to evaluate the PSF (and
                ///derivatives).
                double y,

                ///An array to place the result in.
                ITERABLE result,

                ///Should first order derivatives w.r.t S, D, K be calculated.
                bool calculate_first_deriv=false,

                ///Should secord order derivatives w.r.t. S, D, K be
                ///calculated
                bool calculate_second_deriv=false
        ) const;

        ///The polynomial coefficients of the taylor expansion around x, y.
        double poly_coef(double x, double y, unsigned x_power,
                unsigned y_power) const;

        ///Sets the background for this PSF to the given value
        void set_background(double background) {__bg=background;}

        ///\brief Calculates the integral of the PSF over a rectangle (or its
        ///overlap with a circle).
        ///
        ///Uses the local polynomial approximation around the center of the
        ///rectangle. Subdiving into smaller rectangles if necessary in order
        ///to achive the desired precision.
        ///
        ///The rectangle is defined by center_x - dx/2 < x < center_x + dx/2
        ///and center_y - dy/2 < y < center_y + dy/2. The circle is always
        ///centered on (0, 0) and has the given radius. If the circle radius
        ///is zero, the integral over the full rectangle is calculated.
        double integrate(
            double center_x,
            double center_y,
            double dx,
            double dy,
            double circle_radius=0
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            ,
            bool reset_piece_id=false,
            bool skip_piece=false
#endif
#endif
        ) const
        {
            return integrate(
                center_x,
                center_y,
                dx,
                dy,
                circle_radius,
                false,
                false
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                ,
                reset_piece_id,
                skip_piece
#endif
#endif
            );
        }

        ///\brief Calculates the integral of the PSF over a rectangle (or its
        ///overlap with a circle) with optional derivatives.
        double integrate(double center_x, double center_y,
                double dx, double dy, double circle_radius,
                bool calculate_first_deriv,
                bool calculate_second_deriv
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                ,
                bool reset_piece_id=false, bool skip_piece=false
#endif
#endif
                ) const;

        double last_integrated(SDKDerivative deriv=NO_DERIV) const
        {return __integral_values[deriv];}

        ///Returns the value of S
        double s() const {return 0.5*(__spd+__smd);}

        ///Returns the value of D
        double d() const {return 0.5*(__spd-__smd);}

        ///Returns the value of K
        double k() const {return __k;}

        ///How far left does the PSF extend.
        double min_x() const {return -Core::Inf;}

        ///How far right does the PSF extend.
        double max_x() const {return Core::Inf;}

        ///How far down does the PSF extend.
        double min_y() const {return -Core::Inf;}

        ///How far up does the PSF extend.
        double max_y() const {return Core::Inf;}

        ///Virtual destructor since virtual functions present.
        virtual ~EllipticalGaussian() {}
    };

    template<class ITERABLE>
    void EllipticalGaussian::evaluate(double x, double y,
                                         ITERABLE result,
                                         bool calculate_first_deriv,
                                         bool calculate_second_deriv) const
    {
        result[NO_DERIV] = exp(exp_argument(x,y));
        if(calculate_first_deriv || calculate_second_deriv) {
            double dS_coef = x * x + y * y,
                   dD_coef = x * x - y * y,
                   dK_coef = x * y;
            if(calculate_first_deriv) {
                result[S_DERIV] = dS_coef * result[NO_DERIV];
                result[D_DERIV] = dD_coef * result[NO_DERIV];
                result[K_DERIV] = dK_coef * result[NO_DERIV];
            }
            if(calculate_second_deriv) {
                result[SS_DERIV] = dS_coef * dS_coef * result[NO_DERIV];
                result[SD_DERIV] = dS_coef * dD_coef * result[NO_DERIV];
                result[SK_DERIV] = dS_coef * dK_coef * result[NO_DERIV];

                result[DD_DERIV] = dD_coef * dD_coef * result[NO_DERIV];
                result[DK_DERIV] = dD_coef * dK_coef * result[NO_DERIV];

                result[KK_DERIV] = dK_coef * dK_coef * result[NO_DERIV];
            }
        }
        result[0] += __bg;
        return;
    }

} //End PSF namespace.

#endif
