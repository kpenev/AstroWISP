/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of
 * FitPSF.
 *
 * \ingroup FitPSF
 */

#include "../Core/CInterface.h"

extern "C" {
    ///Opaque struct to cast to/from FitPSF::Config.
    struct LIB_PUBLIC FittingConfiguration;

    ///Create an object for holding the configuration for PSF fitting.
    LIB_PUBLIC FittingConfiguration *create_psffit_configuration(
    );

    ///\brief Destroy a configuration previously created by
    ///create_psffit_configuration()
    LIB_PUBLIC void destroy_psffit_configuration(
        ///The configuration to destroy.
        FittingConfiguration *configuration
    );

    ///Update the configuration for PSF fitting.
    LIB_PUBLIC void update_psffit_configuration(
        ///The configuration to update.
        FittingConfiguration *target_configuration,

        ///Are we deing PRF fittindg (in liu of PSF fitting)?
        bool prffit,

        ///Alternating <parameter name>, <parameter value> pairs, with both
        ///etries being of type char* type.
        ...
    );

};//End extern "C"
