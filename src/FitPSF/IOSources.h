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

#include <vector>
#include <istream>
#include <string>
#include <sstream>
#include <cctype>

namespace FitPSF {
    ///Convenient interface to a single section of FitPSF input source lists.
    class LIB_LOCAL IOSources {
        private:
            ///The locations of the sources for PSF fitting.
            std::vector<Core::SourceLocation> __locations;

            ///See psf_terms() method.
            Eigen::MatrixXd __psf_terms;

            ///See enabled() method.
            std::vector<bool> __enabled;

            std::string
                ///The FITS filename to use for PSF fitting.
                __fits_fname,

                ///The name of the output file to use for the given sources.
                __output_fname,

                ///The name of the file to save source assignemnt info to.
                __source_assignment_fname;

            ///Is this the last section.
            bool __last;

        public:
            ///\brief Construct from an array of sources.
            IOSources(
                ///The FITS filename these sources are contained in.
                const char *fits_fname,

                ///The IDs to assign to these sources.
                char **source_ids,

                ///\brief The coordinates of the sources within the image.
                ///
                ///Should have 2*num_sources entries x0, y0, x1, y1, ...
                const double *source_coordinates,

                ///The terms in the expansion of the PSF map.
                ///
                ///Information about the sources organized in equal sized
                ///columns. The first num_terms entries are the values of
                ///the expansion terms for the first source, followed by the
                ///expansion terms for the second source, etc.
                const double *psf_terms,

                ///Boolean flag for each source indicating whether the source is
                ///allowed to participate in fitting for the PSF shape. All
                ///sources participate in flux fitting.
                const bool *enabled,

                ///How many sources are in column_data.
                unsigned long num_sources,

                ///How many terms are in the PSF expansion.
                unsigned long num_terms
            );

            ///The locations of the sources for PSF fitting.
            const std::vector<Core::SourceLocation> &locations() const
            {return __locations;}

            ///Terms the PSF parameters are linear functions of.
            ///
            ///The terms for source i are the i-th column vector of the returned
            ///matrix.
            const Eigen::MatrixXd &psf_terms() const
            {return __psf_terms;}

            ///Is each source allowed to participate in PSF shape fitting?
            const std::vector<bool> &enabled() const
            {return __enabled;}

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
