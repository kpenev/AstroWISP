/**\file
 *
 * \brief Declares Input/Output interface for the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#ifndef __FIT_PSF_IO_H
#define __FIT_PSF_IO_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../PSF/Typedefs.h"
#include "../Core/SourceLocation.h"
#include "../Core/Typedefs.h"
#include "../Core/Error.h"
#include "../Core/NaN.h"

#include <list>
#include <istream>
#include <string>
#include <sstream>
#include <cctype>

namespace FitPSF {
    ///\brief Parser for a single section of FitPSF input source lists.
    ///
    ///The input file should consist of sections formatted like:
    ///[\<FITS filename\> \<output filename\>]
    ///colvalue1 colvalue2 ...
    ///
    ///At least three columns must be included: ID, x and y.
    ///
    ///White space is ignored (and hence not allowed in filenames).
    class LIB_LOCAL IOSources {
        private:
            ///The locations of the sources for PSF fitting.
            std::list<Core::SourceLocation> __locations;

            ///All non source-ID columns in the same order as the sources.
            PSF::MapVarListType __columns;

            std::string
                ///The FITS filename to use for PSF fitting.
                __fits_fname,

                ///The name of the output file to use for the given sources.
                __output_fname,

                ///The name of the file to save source assignemnt info to.
                __source_assignment_fname;

            ///Is this the last section.
            bool __last;

            ///Check for correct format and read FITS and output filenames.
            void read_filenames(
                ///The stream to read from.
                std::istream &input_stream
            );

            ///Read the column data, initializing __locations along the way.
            void read_column_data(
                ///The stream to read from.
                std::istream &input_stream,

                ///The names of the input columns.
                const std::list<std::string> &column_names,

                ///The location to fill with the data (each list contains
                ///the values for a column).
                /// - Should contain the correct number of columns on
                ///   input.
                /// - Data is appended to the end of the columns so they
                ///   should probably be empty.
                /// - The ID column should not have an entry.
                std::vector< std::list<double> > &column_data
            );

            ///\brief Use the contents of __columns to set the (x, y)
            ///coordinates of  __locations.
            void set_source_coordinates();

        public:
            ///\brief Read into self a single section (one FITS frames's
            ///sources) from the given stream.
            IOSources(
                ///The stream to read from.
                std::istream &input_stream,

                ///The names of the columns in each section.
                const std::list<std::string> &column_names
            );

            ///\brief Construct from an array of sources.
            IOSources(
                ///The FITS filename these sources are contained in.
                const char *fits_fname,

                ///The IDs to assign to these sources.
                char **source_ids,

                ///The information about the sources organized in equal sized
                ///columns. The meaning and the order of the columns is
                ///specified by column_names. The first num_sources entries are
                ///the values of first column, followed by the second column
                ///etc.
                const double *column_data,

                ///The names of the columns.
                char **column_names,

                ///How many sources are in column_data.
                unsigned long num_sources,

                ///How many columns.
                unsigned long num_columns
            );

            ///The locations of the sources for PSF fitting.
            const std::list<Core::SourceLocation> &locations() const
            {return __locations;}

            ///All non source-ID columns in the same order as the sources.
            const PSF::MapVarListType &columns() const
            {return __columns;}

            ///The FITS filename to use for PSF fitting.
            const std::string &fits_fname() const
            {return __fits_fname;}

            ///The output filename to save the PSF fitting results in.
            const std::string &output_fname() const
            {return __output_fname;}

            ///The name of the file to save source assignemnt info to.
            const std::string &source_assignment_fname() const
            {return __source_assignment_fname;}

            ///Was this the last section in the input stream?
            bool last() const {return __last;}
    };
}

#endif
