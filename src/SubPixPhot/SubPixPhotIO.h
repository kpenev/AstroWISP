/**\file
 *
 * \brief Defines I/O utilities for the subpixphot executable.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SUBPIX_PHOT_IO_H
#define __SUBPIX_PHOT_IO_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SourceInput.h"
#include "../Core/Source.h"
#include "../Core/PhotColumns.h"
#include <iostream>
#include <valarray>
#include <list>
#include <iomanip>
#include <cmath>
#include <fstream>

namespace SubPixPhot {

    ///Outputs the sources to stdout
    LIB_LOCAL void output_to_stdout(
        ///The sources to output.
        const std::list<IO::OutputSDKSource>&   sources,

        ///The columns to output.
        const std::list<Phot::Columns>&         columns,

        ///The magnitude that corresponds to a flux of 1ADU
        double                                  mag_1ADU,

        ///The gain with which the image was taken.
        double                                  gain,

        ///The stream to write the output to.
        std::ostream&                           os = std::cout
    );

    ///\brief Output a header to the given stream with version and command 
    ///line information.
    ///
    ///Includes the version of HATpipe if HATP_SVN is defined, the version of
    ///SubPixPhot if the version_string argument is not an empty string, and 
    ///the  command line with which the command was invoked if argc is 
    ///non-zero, in which case argc and argv should be the same arguments 
    ///that main() received.
    LIB_LOCAL void write_header(
        ///The stream to output the header to.
        std::ostream      &os,

        ///The number of command line arguments used.
        int               argc = 0,

        ///The actual command line arguments.
        char**            argv = NULL,

        ///A string which identifies the version of the tool.
        const std::string &version_string = ""
    );

    ///\brief Reads the next source from source_input into the last element 
    ///of sources performing all necessary checks.
    template<class SOURCE_TYPE>
        static void checked_read_last_source(
            SourceInput &source_input,
            std::list<SOURCE_TYPE> &sources
        )
        {
            source_input >> sources.back();
            if(!source_input) { 
                sources.pop_back();
                if(source_input.fail()) {
                    std::ostringstream msg;
                    msg << "Bad source format found for source number "
                        << sources.size();
                    throw Error::IO(msg.str());
                }
            }
        }

    ///\brief Fills the source list with the sources contained in the given 
    ///stream initialized suitably for photometry with the given number of 
    ///apertures. 
    template<class SOURCE_TYPE>
        LIB_PUBLIC void read_sources(
            SourceInput &source_input,
            std::list<SOURCE_TYPE> &sources, 
            unsigned num_apertures
        )
        {
            while(!source_input.eof()) {
#ifdef TRACK_PROGRESS
                std::cerr << "Reading source." << std::endl;
#endif
                sources.push_back(SOURCE_TYPE());
                sources.back().set_num_apertures(num_apertures);
                checked_read_last_source(source_input, sources); 
#ifdef TRACK_PROGRESS
                std::cerr << "Read " << sources.size() << " sources."
                          << std::endl;
#endif
            }
        }

    ///\brief Fills the source list with the sources contained in the given 
    ///stream, SOURCE_TYPE must have a default constructor.
    template<class SOURCE_TYPE>
        LIB_PUBLIC void read_sources(SourceInput &source_input,
                                     std::list<SOURCE_TYPE> &sources)
        {
            while(!source_input.eof()) {
#ifdef TRACK_PROGRESS
                std::cerr << "Reading source." << std::endl;
#endif
                sources.push_back(SOURCE_TYPE());
                checked_read_last_source(source_input, sources); 
#ifdef TRACK_PROGRESS
                std::cerr << "Read " << sources.size() << " sources." << std::endl;
#endif
            }
        }

    ///Reads in polynomial expansion coefficients of S,D, K from a file.
    std::list<double> LIB_PUBLIC read_sdk_coef(
        ///The name of the file to read from.
        const std::string &fname
    );

} //End SubPixPhot namespace.

#endif
