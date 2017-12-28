/**\file
 *
 * \brief Declare a class describing the splitting of areas into smaller
 * pieces.
 *
 * \ingroup PSF
 */

#include "../Core/SharedLibraryExportMacros.h"

namespace PSF {

    ///\brief Describes a splitting of a rectangular area into smaller 
    ///pieces.
    ///
    ///\ingroup PSF
    class LIB_LOCAL Split {
    private:
        int __x_split,///< The number of pieces to split the area in along x.
            __y_split;///< The number of pieces to split the area in along y.

        double __sub_x0, ///< The x coordinate center of the first piece.
               __sub_y0, ///< The y coordinate center of the first piece.
               __sub_dx, ///< The x size of the pieces.
               __sub_dy; ///< The y size of the pieces.

        bool __split; ///< Is any splitting necessary.
    public:
        Split(double spd,///< S+D
                double smd,///< S-D
                double k, ///< K
                double x0, ///< The x coondinate of the center of the area.
                double y0, ///< The y coondinate of the center of the area.
                double dx, ///< The x size of the pieces.
                double dy, ///< The y size of the pieces.
                
                ///The maximum value allowed for any exponent argument.
                double max_exp_coef);

        ///Is any splitting necessary.
        operator bool() {return __split;}

        ///The number of pieces to split the area in along x.
        int x_split() {return __x_split;}

        ///The number of pieces to split the area in along y.
        int y_split() {return __y_split;}

        ///The number of pieces the region was splin into.
        int num_pieces() {return __x_split*__y_split;}

        ///The x coordinate center of the first piece.
        double sub_x0() {return __sub_x0;}

        ///The y coordinate center of the first piece.
        double sub_y0() {return __sub_y0;}

        ///The x size of the pieces.
        double sub_dx() {return __sub_dx;}

        ///The y size of the pieces.
        double sub_dy() {return __sub_dy;} 
    };

} //End PSF namespace.
