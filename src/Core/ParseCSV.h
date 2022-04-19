/**\file
 * \brief Declarations of functions for parsing strings of comma separated
 * values to lists.
 *
 * \ingroup IO
 */

#ifndef __PARSE_CSV_H
#define __PARSE_CSV_H

#include "../Core/SharedLibraryExportMacros.h"
#include "PhotColumns.h"
#include "Typedefs.h"
#include "Error.h"
#include <list>
#include <string>

namespace Core {

    ///\brief Pareses a string of comma separated integer values into a list.
    ///
    ///See parse_real_list() for a description of the arguments.
    LIB_LOCAL std::list<int> parse_int_list(
        const std::string &csv,
        const std::string &optname,
        unsigned min_count,
        unsigned max_count
    );

    ///Parses a string of comma separated real values into a list.
    LIB_LOCAL RealList parse_real_list(
        ///The list to parse
        const std::string &csv,

        ///The name of the option (only used for error message)
        const std::string &optname,

        ///Minimum required number of values (exception thrown if not met).
        unsigned min_count,

        ///Maximum required number of values (exception thrown if exceeded).
        unsigned max_count
    );

    ///Parses a string of comma separated sub-strings into a list.
    LIB_LOCAL StringList parse_string_list(
        ///The list to parse.
        const std::string &csv,
        ///The name of the option (only used for error message)
        const std::string &optname,

        ///Minimum required number of values (exception thrown if not met).
        unsigned min_count,

        ///Maximum required number of values (exception thrown if exceeded).
        unsigned max_count
    );

    ///\brief Parses a string of comma separated column names into a list,
    ///replicating per-aperture columns appropriately.
    LIB_LOCAL ColumnList parse_column_list(

        ///The string to parse.
        const std::string &csv,

        ///How many copies of per-aperture columns to make.
        unsigned num_apertures,

        ///The name of the option being parsed (only used for error message).
        const std::string &optname,

        ///Instead of throwing an exception assign Phot::unknown to any
        ///unrecognized column name.
        bool allow_unknown=true
    );

    ///\brief Perses a string of comma separated values into a list.
    ///
    ///See parse_real_list() for a description of the arguments.
    template<class VAL_TYPE>
        LIB_LOCAL void parse_csv_list(const std::string &csv,
                                      const std::string &optname,
                                      unsigned min_count,
                                      unsigned max_count,
                                      std::list<VAL_TYPE> &result)
        {
            std::istringstream csv_stream(csv);
            VAL_TYPE value;
            while(csv_stream) {
                csv_stream >> value;
                if(csv_stream) result.push_back(value);
                if(
                    csv_stream.eof()
                    &&
                    result.size()>=min_count
                    &&
                    result.size()<=max_count
                )
                    return;
                if(csv_stream.get() != ',') break;
            }
            std::ostringstream error_msg;
            error_msg << "Malformatted "
                << optname
                << " option: "
                << csv
                << "expected comma separated list of "
                << "at least " << min_count << " "
                << "and "
                << "at most " << max_count << " "
                << "real values.";
            throw Error::CommandLine(error_msg.str());
        }

} //End IO namespace.

#endif
