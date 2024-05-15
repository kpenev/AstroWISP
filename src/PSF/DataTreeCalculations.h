/**\file
 *
 * \brief Declare functions that calculate some entries of I/O data trees
 * from others.
 *
 * \ingroup IO
 */

#ifndef __DATA_TREE_CALCULATIONS_H
#define __DATA_TREE_CALCULATIONS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../IO/H5IODataTree.h"
#include "../IO/OutputArray.h"
#include "../PSF/EllipticalGaussianMap.h"
#include "../PSF/MapSourceContainer.h"

#include <cmath>

namespace PSF {

    ///\brief Fill the PSF fluxes in a data tree (computed if necessary).
    ///
    ///If the fluxes are not directly available, magnitudes and magnitude
    ///zero points must be.
    LIB_PUBLIC void fill_psf_fluxes(
        IO::H5IODataTree &data,
        const std::string &data_tree_image_id=""
    );

    ///\brief Fill The PSF amplitudes in a data tree (computed if necessary).
    ///
    ///If the amplitudes are not already available, both the PSF shape
    ///information should be and some form of flux information: fluxes or
    ///magnitudes + magnitude zero point. If fluxes are not available (as
    ///well as amplitudes) they will also be filled.
    LIB_PUBLIC void fill_psf_amplitudes(
        IO::H5IODataTree &data,
        const std::string &data_tree_image_id=""
    );

} //End IO namespace.

#endif
