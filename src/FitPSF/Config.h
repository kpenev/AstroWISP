/**\file
 *
 * \brief Declares the command line parser for the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#ifndef __CONFIG_H
#define __CONFIG_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../PSF/Grid.h"
#include "../PSF/Typedefs.h"
#include "../Background/Annulus.h"
#include "../IO/CommandLineConfig.h"
#include "../PSF/CommandLineUtil.h"
#include "../Background/CommandLineUtil.h"
#include "../Core/CommandLineUtil.h"
#include "../Core/Error.h"
#include "../Core/Typedefs.h"

namespace FitPSF {

    ///\brief Default configuration from file but overwritten by command line
    ///options.
    ///
    ///\ingroup FitPSF
    class LIB_PUBLIC Config : public IO::CommandLineConfig {
    private:
        ///Describes the available command line options.
        void describe_options();

        ///\brief The part of the help describing the usage and purpose (no
        ///options).
        std::string usage_help(
            ///The name of the executable invoked.
            const std::string &prog_name
        ) const;

    public:
        ///\brief Checks for consistency between the command line options and
        ///perform some finalizing configuration.
        ///
        ///Throws an exception if some inconsistency is detected.
        void check_and_finalize();

        ///Parse the command line.
        Config(
            ///The number of arguments on the command line
            ///(+1 for the executable)
            int argc,

            ///A C style array of the actual command line arguments.
            char **argv
        )
        {if(argc > 0) {parse(argc, argv); check_and_finalize();}}
    }; //End Config class.

} //End FitPSF namespace.

#endif
