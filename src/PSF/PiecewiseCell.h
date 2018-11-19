/**\file
 * \brief Defines a base class for the cells comprising a Piecewise.
 */

#ifndef __PIECEWISE_PSF_CELL_H
#define __PIECEWISE_PSF_CELL_H

#include "../Core/SharedLibraryExportMacros.h"
#include <cmath>
#include <valarray>

namespace PSF {

    ///\brief An abstract base class for all cell of a Piecesiwe.
    class LIB_LOCAL PiecewiseCell {
    private:
        ///The horizontal size of the cell.
        double __horizontal_size,

               ///The vertical size of the cell.
               __vertical_size;
    public:
        ///Create a cell with the given size.
        PiecewiseCell(
            ///The horizental size of the new cell.
            double horizontal_size,

            ///The vertical size of the new cell.
            double vertical_size
        ) :
            __horizontal_size(horizontal_size), 
            __vertical_size(vertical_size)
        {}

        ///\brief Calculate the integral of the intensity over a rectangle 
        ///that is contained entirely within the cell.
        ///
        /// \image html PiecewisePSFCell_integrate_rectangle_dddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_rectangle(
            ///The left boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to the left cell
            ///boundary.
            double xmin,

            ///The right boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to the left cell
            ///boundary.
            double xmax,

            ///The bottom boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double ymin,

            ///The top boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double ymax
        ) const =0;

        ///\brief Calculate the integral of the intensity over a rectangle 
        ///that is contained entirely within the cell.
        ///
        /// \image html PiecewisePSFCell_integrate_rectangle_dddd.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_rectangle(
            ///The left boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to the left cell
            ///boundary.
            double xmin,

            ///The right boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to the left cell
            ///boundary.
            double xmax,

            ///The bottom boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double ymin,

            ///The top boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double ymax,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const =0;

        ///\brief Integrate the intensity over a region bounded by one 
        ///vertical, two horizontal lines and an arc that intersects the two 
        ///horizontal lines.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_dddddd.png "The area being integrated along with the meaning of all function arguments.".
        virtual double integrate_hcircle_piece(
            ///The bottom boundary of the region.
            double ymin,

            ///The top boundary of the region.
            double ymax,

            ///The straight vertical boundary of the region.
            double xbound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const = 0;

        ///\brief Integrate the intensity over a region bounded by one
        ///horizontal, two vertical lines and an arc that intersects the two
        ///horizontal lines.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_dddddd.png "The area being integrated along with the meaning of all function arguments.".
        virtual double integrate_vcircle_piece(
            ///The left boundary of the region.
            double xmin,

            ///The right boundary of the region.
            double xmax,

            ///The straight horizontal boundary of the region.
            double ybound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const = 0;

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell partially in the vertical direction.
        ///
        /// \image html PiecewisePSFCell_integrate_partial_vspan_dddb.png "The area being integrated along with the meaning of all function arguments. Left: up=true, right: up=false".
        virtual double integrate_partial_vspan(
            ///The left boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_min,

            ///The right boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_max,

            ///The bottom/top y boundary of the rectangle if top is
            ///true/false. The other boundary is the top/bottom cell
            ///boundary.  Defined relative to the bottom cell boundary.
            double y,

            ///Wether the rectangle extends upward of y.
            bool up
        ) const
        {return integrate_rectangle(x_min, x_max, (up ? y : 0),
                                    (up ? __vertical_size : y));}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell partially in the vertical direction.
        ///
        /// \image html PiecewisePSFCell_integrate_partial_vspan_dddb.png "The area being integrated along with the meaning of all function arguments. Left: up=true, right: up=false".
        virtual std::valarray<double> integrate_partial_vspan(
            ///The left boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_min,

            ///The right boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_max,

            ///The bottom/top y boundary of the rectangle if top is
            ///true/false. The other boundary is the top/bottom cell
            ///boundary.  Defined relative to the bottom cell boundary.
            double y,

            ///Wether the rectangle extends upward of y.
            bool up,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_rectangle(x_min, x_max, (up ? y : 0),
                                    (up ? __vertical_size : y), coef_sets);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell partially in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_partial_hspan_dddb.png "The area being integrated along with the meaning of all function arguments. Left: right=false, right: right=true".
        virtual double integrate_partial_hspan(
            ///The bottom boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y_min,

            ///The top boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y_max,

            ///The left/right y boundary of the rectangle if right is
            ///true/false. The other boundary is the right/left cell
            ///boundary.  Defined relative to the left cell boundary.
            double x,

            ///Wether the rectangle extends rightward of x.
            bool right
        ) const
        {return integrate_rectangle((right ? x : 0),
                                    (right ? __horizontal_size : x),
                                    y_min,
                                    y_max);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell partially in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_partial_hspan_dddb.png "The area being integrated along with the meaning of all function arguments. Left: right=false, right: right=true".
        virtual std::valarray<double> integrate_partial_hspan(
            ///The bottom boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y_min,

            ///The top boundary of the rectangle. Must lie between the 
            ///vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y_max,

            ///The left/right y boundary of the rectangle if right is
            ///true/false. The other boundary is the right/left cell
            ///boundary.  Defined relative to the left cell boundary.
            double x,

            ///Wether the rectangle extends rightward of x.
            bool right,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_rectangle((right ? x : 0),
                                    (right ? __horizontal_size : x),
                                    y_min,
                                    y_max,
                                    coef_sets);}

        ///Calculate the integral of the intensity over the enitire cell.
        virtual double integrate() const
        {
            return integrate_rectangle(0,
                                       __horizontal_size,
                                       0,
                                       __vertical_size);
        }

        ///Calculate the integral of the intensity over the enitire cell.
        virtual std::valarray<double> integrate(
            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets) const
        {
            return integrate_rectangle(0,
                                       __horizontal_size,
                                       0,
                                       __vertical_size,
                                       coef_sets);
        }

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the vertical direction..
        ///
        /// \image html PiecewisePSFCell_integrate_vspan_dd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vspan(
                ///The left boundary of the rectangle. Must lie between the 
                ///horizontal cell boundaries. Defined relative to 
                ///the left cell boundary.
                double x_min,

                ///The right boundary of the rectangle. Must lie between the 
                ///horizontal cell boundaries. Defined relative to 
                ///the left cell boundary.
                double x_max) const
        {return integrate_partial_vspan(x_min, x_max, 0, true);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the vertical direction..
        ///
        /// \image html PiecewisePSFCell_integrate_vspan_dd.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_vspan(
            ///The left boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_min,

            ///The right boundary of the rectangle. Must lie between the 
            ///horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x_max,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_partial_vspan(x_min, x_max, 0, true, coef_sets);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_hspan_dd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hspan(
            ///The bottom boundary of the rectangle. Must lie between the
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double y_min,

            ///The top boundary of the rectangle. Must lie between the
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double y_max
        ) const
        {return integrate_partial_hspan(y_min, y_max, 0, true);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_hspan_dd.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_hspan(
            ///The bottom boundary of the rectangle. Must lie between the
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double y_min,

            ///The top boundary of the rectangle. Must lie between the
            ///vertical cell boundaries. Defined relative to the bottom cell
            ///boundary.
            double y_max,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_partial_hspan(y_min, y_max, 0, true, coef_sets);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the vertical direction..
        ///
        /// \image html PiecewisePSFCell_integrate_vspan_db.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vspan(
            ///The left/right boundary of the rectangle if right is true/false
            ///the other boundary is the right/left cell boundary. Must lie
            ///between the horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x,

            ///Does the rectangle extend rightward of the x argument?
            bool right
        ) const
        {return integrate_vspan((right ? x : 0),
                (right ? __horizontal_size : x));}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the vertical direction..
        ///
        /// \image html PiecewisePSFCell_integrate_vspan_db.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_vspan(
            ///The left/right boundary of the rectangle if right is true/false
            ///the other boundary is the right/left cell boundary. Must lie
            ///between the horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x,

            ///Does the rectangle extend rightward of the x argument?
            bool right,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_vspan((right ? x : 0),
                (right ? __horizontal_size : x), coef_sets);}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_hspan_db.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hspan(
            ///The top/bottom boundary of the rectangle if top is false/true
            ///the other boundary is the bottom/top cell boundary. Must lie
            ///between the vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y,

            ///Does the rectangle extend upward of the y argument?
            bool up
        ) const
        {return integrate_hspan((up ? y : 0), (up ? __vertical_size : y));}

        ///\brief Calculate the integral of the intensity over a rectangle
        ///spanning the cell fully in the horizontal direction..
        ///
        /// \image html PiecewisePSFCell_integrate_hspan_db.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_hspan(
            ///The top/bottom boundary of the rectangle if top is false/true
            ///the other boundary is the bottom/top cell boundary. Must lie
            ///between the vertical cell boundaries. Defined relative to 
            ///the bottom cell boundary.
            double y,

            ///Does the rectangle extend upward of the y argument?
            bool up,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_hspan((up ? y : 0), (up ? __vertical_size : y),
                coef_sets);}

        ///\brief Calculate the integral of the intensity over a rectangle.
        ///
        /// \image html PiecewisePSFCell_integrate_rectangle_ddbb.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_rectangle(
            ///The left/right boundary of the rectangle if right is true/false
            ///the other boundary is the right/left cell boundary. Must lie
            ///between the horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x,

            ///The top/bottom boundary of the rectangle if top is false/true
            ///the other boundary is the bottom/top cell boundary. Must lie
            ///between the vertical cell boundaries. Defined relative to the
            ///bottom cell boundary.
            double y,

            ///Does the rectangle extend rightward of the x argument?
            bool right,

            ///Does the rectangle extend upward of the y argument?		
            bool up
        ) const
        {return integrate_rectangle(
                (right ? x : 0), (right ? __horizontal_size : x),
                (up ? y : 0), (up ? __vertical_size : y));}

        ///\brief Calculate the integral of the intensity over a rectangle.
        ///
        /// \image html PiecewisePSFCell_integrate_rectangle_ddbb.png "The area being integrated along with the meaning of all function arguments."
        virtual std::valarray<double> integrate_rectangle(
            ///The left/right boundary of the rectangle if right is true/false
            ///the other boundary is the right/left cell boundary. Must lie
            ///between the horizontal cell boundaries. Defined relative to 
            ///the left cell boundary.
            double x,

            ///The top/bottom boundary of the rectangle if top is false/true
            ///the other boundary is the bottom/top cell boundary. Must lie
            ///between the vertical cell boundaries. Defined relative to the
            ///bottom cell boundary.
            double y,

            ///Does the rectangle extend rightward of the x argument?
            bool right,

            ///Does the rectangle extend upward of the y argument?		
            bool up,

            ///The sets of coefficients to compute the integral for. Each set
            ///must contain 16 coefficients and the sets are concatenated.
            const std::valarray<double> &coef_sets
        ) const
        {return integrate_rectangle(
                (right ? x : 0), (right ? __horizontal_size : x),
                (up ? y : 0), (up ? __vertical_size : y), coef_sets);}

        ///\brief Same as integrate_hcircle_piece(double, double, double, 
        ///double, double, double) but with either the top or the bottom 
        ///boundary coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_dddddb.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hcircle_piece(
            ///The top/bottom boundary of the region if up is false/true.
            double ybound,

            ///The straight vertical boundary of the region.
            double xbound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y,

            ///Wether the region extends upwards from ybound.
            bool up
        ) const
        {return integrate_hcircle_piece( (up ? ybound : 0),
                (up ? __vertical_size : ybound), xbound, radius,
                circle_x, circle_y);}

        ///\brief Same as integrate_hcircle_piece(double, double, double, 
        ///double, double, double) but with the vertical boundary coinciding 
        ///with the cell wall inside the circle of the boundinc arc.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_ddddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hcircle_piece(
            ///The bottom boundary of the region.
            double ymin,

            ///The top boundary of the region.
            double ymax,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            return integrate_hcircle_piece(
                ymin,
                ymax,
                (circle_x<=0 ? 0 : __horizontal_size),
                radius,
                circle_x,
                circle_y
            );
        }

        ///\brief Same as integrate_hcircle_piece(double, double, double, 
        ///double, double, double) but with both the top and the bottom 
        ///boundary coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_dddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hcircle_piece(
            ///The straight vertical boundary of the region.
            double xbound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            return integrate_hcircle_piece(0,
                                           xbound,
                                           radius,
                                           circle_x,
                                           circle_y,
                                           true);
        }

        ///\brief Same as integrate_hcircle_piece(double, double, double, 
        ///double, double, double) but with all three straight boundaries 
        ///coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_ddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_hcircle_piece(
            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            return integrate_hcircle_piece(
                (circle_x<=0 ? 0 : __horizontal_size),
                radius,
                circle_x,
                circle_y
            );
        }

        ///\brief Same as integrate_vcircle_piece(double, double, double, 
        ///double, double, double) but with either the left or the right 
        ///boundary coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_dddddb.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vcircle_piece(
            ///The left/right boundary of the region if right is true/false.
            double xbound,

            ///The straight horizontal boundary of the region.
            double ybound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y,

            ///Wether the region extends rightward from xbound.
            bool right
        ) const
        {
            return integrate_vcircle_piece(
                (right ? xbound : 0),
                (right ? __horizontal_size : xbound),
                ybound,
                radius,
                circle_x,
                circle_y
            );
        }

        ///\brief Same as integrate_vcircle_piece(double, double, double, 
        ///double, double, double) but with the horizontal boundary 
        ///coinciding with the cell wall inside the circle of the boundinc 
        ///arc.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_ddddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vcircle_piece(
            ///The bottom boundary of the region.
            double xmin,

            ///The top boundary of the region.
            double xmax,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y) const
        {
            return integrate_vcircle_piece(
                xmin,
                xmax,
                (circle_y<=0 ? 0 : __vertical_size),
                radius,
                circle_x,
                circle_y
            );
        }

        ///\brief Same as integrate_vcircle_piece(double, double, double, 
        ///double, double, double) but with both the left and the right 
        ///boundary coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_dddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vcircle_piece(
            ///The straight horizontal boundary of the region.
            double ybound,

            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            return integrate_vcircle_piece(0,
                                           ybound,
                                           radius,
                                           circle_x,
                                           circle_y,
                                           true);
        }

        ///\brief Same as integrate_vcircle_piece(double, double, double, 
        ///double, double, double) but the horizontal as well as one of the 
        ///vertical boundaries coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_ddddb.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vcircle_piece(
                ///The vertical boundary of the region not coinciding with a cell
                ///wall.
                double xbound,

                ///The radius of the circle the bounding arc is part of.
                double radius,

                ///The x coordinate of the circle center relative to the cell
                ///center.
                double circle_x,

                ///The y coordinate of the circle center relative to the cell
                ///center.
                double circle_y,

                ///Whether the region extends to the right of xbound.
                bool right) const
        {
            return integrate_vcircle_piece(xbound,
                                           (circle_y>0 ? __vertical_size : 0),
                                           radius, circle_x, circle_y,
                                           right);
        }

        ///\brief Same as integrate_vcircle_piece(double, double, double, 
        ///double, double, double) but with all three straight boundaries 
        ///coinciding with the cell's.
        ///
        /// \image html PiecewisePSFCell_integrate_vcircle_piece_ddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_vcircle_piece(
            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            return integrate_vcircle_piece(
                (circle_y<=0 ? 0 : __vertical_size),
                radius,
                circle_x,
                circle_y
            );
        }

        ///\brief Integrates the intensity over a region with boundaries two
        ///adjacent cell walls and an arc of a circle intersecting both of 
        ///them.
        ///
        /// \image html PiecewisePSFCell_integrate_wedge_ddd.png "The area being integrated along with the meaning of all function arguments."
        virtual double integrate_wedge(
            ///The radius of the circle the bounding arc is part of.
            double radius,

            ///The x coordinate of the circle center relative to the cell
            ///center.
            double circle_x,

            ///The y coordinate of the circle center relative to the cell
            ///center.
            double circle_y
        ) const
        {
            if(circle_x<=0) {
                double ycross=std::abs(
                        std::sqrt(std::pow(radius,2)-std::pow(circle_x,2))
                        -
                        std::abs(circle_y));
                return integrate_hcircle_piece((circle_y<=0 ? 0 : ycross),
                        (circle_y<=0 ? ycross : __vertical_size), radius,
                        circle_x, circle_y);
            } else {
                double ycross=std::abs(
                        std::sqrt(
                            std::pow(radius,2)
                            -
                            std::pow(circle_x-__horizontal_size, 2)
                            )
                        -
                        std::abs(circle_y)),
                       ymin,
                       ymax;
                if(circle_y<=0) {
                    ymin=0;
                    ymax=ycross;
                } else {
                    ymin=ycross;
                    ymax=__vertical_size;
                }
                if(ymin>=ymax) return 0.0;
                return integrate_hcircle_piece(ymin,
                                               ymax,
                                               radius,
                                               circle_x,
                                               circle_y);
            }
        }

        ///\brief Returns the value of the intensity at the given position
        ///relative to the bottom left cell corner.
        virtual double operator()(
            ///The horizontal offset from the left cell wall where the
            ///intensity should be returned
            double x,

            ///The vertical offset from the bottom cell wall where the
            ///intensity should be returned
            double y
        ) const = 0;

        ///\brief Returns the value of the intensity at the given position
        ///relative to the bottom left cell corner.
        virtual std::valarray<double> operator()(
            ///The horizontal offset from the left cell wall where the
            ///intensity should be returned
            double x,

            ///The vertical offset from the bottom cell wall where the
            ///intensity should be returned
            double y,

            ///The sets of coefficients to compute the intensity for. Each
            ///set must contain 16 coefficients and the sets are
            ///concatenated.
            const std::valarray<double> &coef_sets
        ) const = 0;

        ///Return a newly allocated copy of this object.
        virtual PiecewiseCell* clone() const = 0;

        ///The horizontal size of the cell.
        double horizontal_size() const {return __horizontal_size;}

        ///The vertical size of the cell.
        double vertical_size() const {return __vertical_size;}

        ///Virtual destructor.
        virtual ~PiecewiseCell() {}
    }; //End PiecewiseCell class.

} //End PSF namespace.

#endif
