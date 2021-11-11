/**\file
 *
 * \brief Defines the base class of all PSF models.
 *
 * \ingroup PSF
 */

#ifndef __PSF_H
#define __PSF_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include <valarray>
#include <cmath>

namespace PSF {

    ///An abstract parent class for all PSF models.
    class LIB_PUBLIC PSF {
    private:
        ///Returns the points at which the line segment betwenn (x1, y) and
        ///(x2, y) intersects the circle centered at (0,0) with a radius r.
        ///The arguments (left/right)_on_circle indicate that the point at
        ///(x(1/2), y) respectively lie exactly on the circle, in which case
        ///they are not included in the list of intersections returned.
        ///It is assumed that x1<x2, y>0, and x1 and x2 have the same sign.
        ///The return values have the following meaning:
        /// * 0 : the line is entirely outside the circle (or is just
        ///       touching it)
        /// * NaN : the line is entirely within the circle
        /// * positive value : a single point of intersection exists at the
        ///                    returned value
        /// * negative value : two intersections exist at +- the returned
        ///                    value
        double line_circle_intersections(
            double x1,
            double x2,
            double y,
            double r,
            bool left_on_circle=false,
            bool right_on_circle=false
        ) const;

        ///Integrates the PSF over the intersection of a rectangle and a cirle
        ///under the assumption that the bottom side of the rectangle lies
        ///outside the circle. See integrate_overlap for description of
        ///command line arguments. Also assumes the rectangle is entirely
        ///contained within a quadrant.
        double integrate_overlap_bottom_out(
            ///See same name argu ment to integrate_overlap()
            double x1,

            ///See same name argu ment to integrate_overlap()
            double y1,

            ///See same name argu ment to integrate_overlap()
            double x2,

            ///See same name argu ment to integrate_overlap()
            double y2,

            ///See same name argu ment to integrate_overlap()
            double rc,

            ///See same name argu ment to integrate_overlap()
            std::valarray<bool> on_circle,

            std::valarray<double> &intersections
        ) const;

        ///Integrates the PSF over the intersection of a rectangle and a cirle
        ///under the assumption that the bottom side of the rectangle lies
        ///inside the circle. See integrate_overlap for description of command
        ///line arguments. Also assumes the rectangle is entirely contained
        ///within a quadrant.
        double integrate_overlap_bottom_in(
            ///See same name argu ment to integrate_overlap()
            double x1,

            ///See same name argu ment to integrate_overlap()
            double y1,

            ///See same name argu ment to integrate_overlap()
            double x2,

            ///See same name argu ment to integrate_overlap()
            double y2,

            ///See same name argu ment to integrate_overlap()
            double rc,

            ///See same name argu ment to integrate_overlap()
            std::valarray<bool> on_circle,

            std::valarray<double> &intersections
        ) const;

        ///\brief Integrates the PSF over the common overlapping area of a
        ///rectangle and a circle.
        double integrate_overlap(
            ///The left boundary of the rectangle.
            double x1,

            ///The bottom boundary of the rectangle.
            double y1,

            ///The right boundary of the rectangle. Must be > x1.
            double x2,

            ///The top boundary of the rectangle.
            double y2,

            ///The radius of the circle. The circle is centered at (0, 0).
            double rc,

            ///Identifies if rectangle corners, in the order (x1, y1), (x2, y1,
            ///(x2, y2), (x1, y2), lie exactly on the circle.
            const std::valarray<bool> &on_circle=std::valarray<bool>(false, 4)
        ) const;

    protected:
        ///\brief Calculates the integral of the PSF over a rectangle.
        virtual double integrate_rectangle(
            ///The x coordinate of the center of the rectangle.
            double center_x,

            ///The y coordinate of the center of the rectangle.
            double center_y,

            ///The full x size of the rectangle.
            double dx,

            ///The full y size of the rectangle.
            double dy
        ) const =0;

        ///Integrates the PSF over the smallest of the four wedges with the
        //following boundaries:
        /// * the line x=x
        /// * the line y=y
        /// * the circle centered at (0, 0) with a radius=radius
        ///If x is 0 the the left vs right wedge is chosen according to left
        ///Same for y0 and bottom.
        virtual double integrate_wedge(double x, double y, double radius,
                bool left=false, bool bottom=false) const =0;
    public:
        ///Set a precision requirement for integrals. May be ignored by actual
        ///PSF.
        virtual void set_precision(
            ///Relative precision
            double,

            ///Absolute precision
            double
        ) const
        {}

        ///Evaluates the PSF at the given position
        virtual double operator()(double x, double y) const =0;

        ///\brief Calculates the integral of the PSF over (a piece of, if
        ///circle_radius!=0) a rectangle interior to a circle.
        ///
        ///The rectangle is defined by:
        ///center_x - dx/2 < x < center_x + dx/2,
        ///center_y - dy/2 < y < center_y + dy/2
        ///and the circle is centered at 0 and has the given radius
        virtual double integrate(
            ///The x coordinate of the center of the rectangle to integrate
            ///over.
            double center_x,

            ///The y coordinate of the center of the rectangle to integrate
            ///over.
            double center_y,

            ///The width of the rectangle.
            double dx,

            ///The height of the rectangle.
            double dy,

            ///The radius of the circle.
            double circle_radius=0
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            ,
            ///Should counting PSF pieces start from scratch.
            bool reset_piece_id=false,

            ///If true, information about the piece being integrated is not
            ///output.
            bool skip_piece=false
#endif
#endif
        ) const;

        ///Virtual destructor.
        virtual ~PSF() {}
    }; //End PSF class.

} //End PSF namespace.

#endif
