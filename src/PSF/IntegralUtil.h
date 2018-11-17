/**\file
 *
 * \brief Useful functions for calculating PSF integrals.
 *
 * \ingroup SubPixPhot
 * \ingroup FitPSF
 * \ingroup FitSubpix
 */

#ifndef __INTEGRAL_UTIL_H
#define __INTEGRAL_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Typedefs.h"
#include <vector>

namespace PSF {

    ///\brief Extend a vector with powers x^n up to max_pow, assuming it
    ///already contains some.
    LIB_LOCAL void fill_powers(
        ///The vector to extend. Must already contain at lesat two entries: 1
        ///and x.
        std::vector<double> &powers,

        ///The largest power to include in the vector.
        Core::vector_size_type max_pow
    );

    ///Initializes a vector of powers.
    LIB_LOCAL void initialize_powers(
            ///The vector to initialize.
            std::vector<double> &powers,

            ///The value whose powers are to be stored.
            double x,

            ///Pre-fill powers up to at least this value.
            Core::vector_size_type max_pow,

            ///Make sure the vector can hold at least this many values before
            ///it needs to expand its storage. 
            Core::vector_size_type initial_storage
    );

} //End PSF namespace.

#endif 
