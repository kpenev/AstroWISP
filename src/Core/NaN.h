/**\file
 *
 * \brief Defines not-a-number and infinity.
 *
 * \ingroup SubPixPhot
 * \ingroup FitSubpix
 * \ingroup FitPSF
 */

#ifndef __NAN_H
#define __NAN_H

#include "../Core/SharedLibraryExportMacros.h"
#include <limits>

namespace Core {

    ///Not-a-number.
    LIB_PUBLIC const double NaN = std::numeric_limits<double>::quiet_NaN();

    ///Positive infinity
    LIB_PUBLIC const double Inf = std::numeric_limits<double>::infinity();

    LIB_LOCAL inline void silence_unused_nan_warning(double =NaN) {};
    LIB_LOCAL inline void silence_unused_inf_warning(double =Inf) {};

} //End Core namespace.

#endif
