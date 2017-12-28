/**\file
 * \brief Declarations of utilities for command line parsing backgrund
 * information.
 *
 * \ingroup Background
 */

#ifndef __BACKGROUND_COMMAND_LINE_UTIL_H
#define __BACKGROUND_COMMAND_LINE_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Annulus.h"
#include "../Core/ParseCSV.h"
#include "../Core/Error.h"
#include <string>
#include <boost/program_options.hpp>
#include <list>
#include <fstream>
#include <sstream>
#include <cmath>

namespace opt = boost::program_options;

namespace Background {

    ///Parse a background annulus option
    LIB_PUBLIC void validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation)
        Annulus*, 

        ///Workaround (see boost documentation)
        int
    );
} //End Background namespace.

#endif
