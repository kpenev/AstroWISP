/**\file
 *
 * \brief Defines a source class which can be passed directly to a PSF map.
 *
 * \ingroup FitPSF
 */

#ifndef __PSF_MAP_SOURCE_H
#define __PSF_MAP_SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../IO/H5IODataTree.h"
#include "../Core/Source.h"
#include "../PSF/Typedefs.h"
#include <vector>

namespace PSF {

    ///\brief Class for evaluating a collection of smooth dependencies for a
    ///collection of sources.
    class LIB_LOCAL MapSource : public Core::Source {
    private:
        ///The values of the terms participating in the PSF map.
        Eigen::VectorXd __expansion_terms;

    public:
        ///Create the source with the given properties.
        MapSource(
                ///See Source::Source
                const Core::SourceID &id,

                ///See Source::Source
                unsigned num_apertures,

                ///See Source::Source
                double x0,

                ///See Source::Source
                double y0,

                ///See Source::Source
                const Background::Source &background
        ) : Source(id, num_apertures, x0, y0, background) {}

        MapSource() : Source(Core::SourceID())
        {throw Error::Runtime("Using default MapSource constructor");}

        ///Modifiable reference to the terms PSF is a function of.
        Eigen::VectorXd &expansion_terms()
        {return __expansion_terms;}

        ///Constant reference to the terms PSF is a function of.
        const Eigen::VectorXd &expansion_terms() const
        {return __expansion_terms;}
    };

} //End PSF namespace.

#endif
