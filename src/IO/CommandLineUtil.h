/**\file
 * \brief Declarations of utilities for command line parsing of IO related
 * options.
 *
 * \ingroup IO
 */

#ifndef __IO_COMMAND_LINE_UTIL_H
#define __IO_COMMAND_LINE_UTIL_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include "../Core/PhotColumns.h"
#include <string>
#include <boost/program_options.hpp>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>

namespace opt = boost::program_options;

namespace IO {

    ///\brief Throws an error if columns is missing an entry required for 
    ///photometry.
    LIB_LOCAL void verify_input_columns(
        ///A list of the input columns for photometry.
        const std::list<Phot::Columns> &columns,

        ///Whether PSF information is required
        bool need_PSF=true
    );

} //End IO namespace.


#endif
