/**\file
 * \brief Declarations of utilities for command line parsing of PSF
 * information.
 *
 * \ingroup PSF
 */

#ifndef __PSF_COMMAND_LINE_UTIL_H
#define __PSF_COMMAND_LINE_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../IO/parse_grid.h"
#include "EllipticalGaussian.h"
#include "Grid.h"
#include "Typedefs.h"
#include <boost/program_options.hpp>

namespace opt = boost::program_options;

namespace PSF {

    ///Parse a string defining a grid over which piecewise PSFs are defined.
    void LIB_LOCAL validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation).
        Grid*,

        ///Workaround (see boost documentation).
        int
    );

    ///Parse a string identifying a PSF model to use.
    void LIB_LOCAL validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation).
        ModelType*,

        ///Workaround (see boost documentation).
        int
    );

} //End PSF namespace.

#endif
