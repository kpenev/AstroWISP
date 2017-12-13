/**\file
 *
 * \brief Declares the command line parser for the FitPSF tool.
 *
 * \ingroup FitPSF
 */

#ifndef __CONFIG_H
#define __CONFIG_H

#include "../PSF/Grid.h"
#include "../PSF/Typedefs.h"
#include "../PSF/TermCalculator.h"
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
    class Config : public IO::CommandLineConfig {
    private:
        ///Describes the available command line options.
        void describe_options();

        ///\brief The part of the help describing the usage and purpose (no 
        ///options).
        std::string usage_help(const std::string &prog_name) const;

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
            parse(argc, argv);
            if(count("help")==0) check_consistency();
            PSF::EllipticalGaussian::set_default_precision(
                operator[]("psf.sdk.rel-int-precision").as<double>(),
                operator[]("psf.sdk.abs-int-precision").as<double>()
            );
            PSF::EllipticalGaussian::set_default_max_exp_coef(
                operator[]("psf.sdk.max-exp-coef").as<double>()
            );
        }
    }; //End Config class.

} //End FitPSF namespace.

#endif
