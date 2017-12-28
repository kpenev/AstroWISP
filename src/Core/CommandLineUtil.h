/**\file
 * \brief Declarations of utilities for command line parsing common to all
 * tools.
 *
 * \ingroup Core
 */

#ifndef __CORE_COMMAND_LINE_UTIL_H
#define __CORE_COMMAND_LINE_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "ParseCSV.h"
#include "Typedefs.h"
#include "Error.h"
#include <boost/program_options.hpp>
#include <limits>

namespace opt = boost::program_options;

//Validators need to be in the same namespace as their class for Boost to
//find them!
namespace Core {

    ///Parse a commas separated list of real values to RealList.
    LIB_PUBLIC void validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation)
        RealList*, 

        ///Workaround (see boost documentation)
        int
    );

    ///Parse a commas separated list of strings to StringList.
    LIB_PUBLIC void validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation)
        StringList*, 

        ///Workaround (see boost documentation)
        int
    );

    ///Parse a list of column names option
    LIB_PUBLIC void validate(
        ///The value to parse into.
        boost::any& value,

        ///The strings passed to this option on the command line.
        const std::vector<std::string>& option_strings,

        ///Workaround (see boost documentation)
        ColumnList*, 

        ///Workaround (see boost documentation)
        int
    );

} //End Core namespace.


#endif
