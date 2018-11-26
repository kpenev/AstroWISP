/**\file
 *
 * \brief Defines a class containing the configuration with which to run
 * subpixphot.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SUB_PIX_PHOT_CONFIG_H
#define __SUB_PIX_PHOT_CONFIG_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../IO/CommandLineConfig.h"
#include "../Background/Annulus.h"
#include "../Background/CommandLineUtil.h"
#include "../Core/PhotColumns.h"
#include "../Core/Error.h"
#include "../Core/Typedefs.h"
#include "../Core/CommandLineUtil.h"
#include <string>
#include <sstream>

namespace SubPixPhot {

    ///\brief Default configuration from file but overwritten by command line 
    ///options.
    class LIB_PUBLIC Config : public IO::CommandLineConfig {
    private:
        ///Describes the available command line options.
        void describe_options();

        ///The part of the help describing the usage and purpose (no options).
        std::string usage_help(
            ///The name with which the executable was invoked on the command
            ///line.
            const std::string &prog_name
        ) const;

        ///Uses io.sources instead of io.psfmap or io.output, if not specified.
        void apply_fallbacks();

        ///\brief Checks for consistency between the command line options.
        ///
        ///Throws an exception if some inconsistency is detected.
        void check_consistency();

    public:
        ///Parse the command line.
        Config(
            ///The number of arguments on the command line
            ///(+1 for the executable)
            int argc,

            ///A C style array of the actual command line arguments.
            char **argv
        )
        {
            if(argc > 0) {
                parse(argc, argv);

                if(count("help")==0) {
                    apply_fallbacks();
                    check_consistency();
                }
                PSF::EllipticalGaussian::set_default_precision(
                    operator[]("psf.sdk.rel-int-precision").as<double>(),
                    operator[]("psf.sdk.abs-int-precision").as<double>()
                );
                PSF::EllipticalGaussian::set_default_max_exp_coef(
                    operator[]("psf.sdk.max-exp-coef").as<double>()
                );
            }
        }
    };

} //End SubPixPhot namespace.

#endif
