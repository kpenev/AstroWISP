/**\file
 *
 * \brief Some useful typedef statements.
 *
 * \ingroup PSF
 */

#ifndef __PSF_TYPEDEFS_H
#define __PSF_TYPEDEFS_H

#include "../Core/SharedLibraryExportMacros.h"
#include <utility>
#include <string>
#include <valarray>
#include <list>

namespace PSF {
    ///The type to use for PSF map variables (quantities PSF depends on).
    typedef std::pair<std::string,
                      std::valarray<double> > MapVariableType;

    ///List of PSF map variables.
    typedef std::list<MapVariableType> MapVarListType;

    ///The allowed PSF/PRF models.
    enum LIB_PUBLIC ModelType {SDK, BICUBIC, ZERO};

}

#endif
