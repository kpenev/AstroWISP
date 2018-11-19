/**\file
 * \brief Defines a PSF cell with intensity given by a bi-cubic polynomial.
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "PiecewiseCell.h"
#include "CirclePieceIntegral.h"
#include "../Core/Typedefs.h"
#include <valarray>
#include <cassert>

namespace PSF {

    ///\brief A class for Piecewise cells over which the intensity is 
    ///given by a bi-cubic polynomial.
    class LIB_PUBLIC PiecewiseBicubicCell : public PiecewiseCell {
    private:
        ///\brief The sets of polynomial coefficients.
        ///
        ///The coefficient in front of \f$x^m y^n\f$ is __coef[m+4*n].
        std::valarray<double> __coef;

    public:
        ///\brief Create a cell with the given size, but without initializing
        ///the polynomial coefficients.
        PiecewiseBicubicCell(
            ///The horizontal size of the cell to create.
            double horizontal_size,

            ///The vertical size of the cell to create.
            double vertical_size
        ) :
            PiecewiseCell(horizontal_size, vertical_size)
        {}

        ///\brief Create a cell with the given horizontal and vertical sizes and
        ///initialize the coefficients.
        ///
        ///See __coef for the order of the coefficients.
        template<class IteratorType>
        PiecewiseBicubicCell(
            ///The horizental size of the cell to create.
            double horizontal_size,

            ///The vertical size of the cell to create.
            double vertical_size,

            ///An iterator to the first polynomial coefficient.
            IteratorType first_coef,

            ///An iterator to one past the last polynomial coefficient.
            IteratorType last_coef
        ) :
            PiecewiseCell(horizontal_size, vertical_size),
            __coef(16)
        {
            for(unsigned i=0; i<16; ++i) {
                assert(!std::isnan(*first_coef));
                __coef[i]=*first_coef; ++first_coef;
            }
            ( void )last_coef;	// unused in release mode
            assert(first_coef==last_coef);
        }

        ///Copy constructor.
        PiecewiseBicubicCell(
            ///The original cell to copy.
            const PiecewiseBicubicCell &orig
        ) :
            PiecewiseCell(orig), __coef(orig.__coef)
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
        ) const
        {return integrate_rectangle(xmin, xmax, ymin, ymax, __coef)[0];}

        ///\brief Calculate the integral of the intensity over a rectangle 
        ///that is contained entirely within the cell for several sets of
        ///coefficients.
        ///
        /// \image html PiecewisePSFCell_integrate_rectangle_dddd.png "The area being integrated along with the meaning of all function arguments."
        std::valarray<double> integrate_rectangle(
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
        ) const;

        ///\brief Integrate the intensity over a region bounded by one 
        ///vertical, two horizontal lines and an arc that intersects the two 
        ///horizontal lines.
        ///
        /// \image html PiecewisePSFCell_integrate_hcircle_piece_dddddd.png  "The area being integrated along with the meaning of all function arguments.".
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
        ) const;

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
        ) const;

        ///\brief Returns the value of the intensity at the given position
        ///relative to the bottom left cell corner.
        virtual double operator()(
            ///The horizontal offset from the left cell wall where the
            ///intensity should be returned
            double x,

            ///The vertical offset from the bottom cell wall where the
            ///intensity should be returned
            double y
        ) const
        {return operator()(x, y, __coef)[0];};

        ///\brief Returns the value of the intensity at the given position
        ///relative to the bottom left cell corner.
        std::valarray<double> operator()(
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
        ) const;

        ///Return a newly allocated copy of this object.
        virtual PiecewiseCell* clone() const
        {return new PiecewiseBicubicCell(*this);}

        using PiecewiseCell::integrate_rectangle;
        using PiecewiseCell::integrate_hcircle_piece;
        using PiecewiseCell::integrate_vcircle_piece;
    }; //End PiecewiseBicubicCell class.

} //End PSF namespace.
