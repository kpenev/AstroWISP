/**\file
 *
 * \brief Declares a container for holding MapSource objects.
 *
 * \ingroup PSF
 */

#ifndef __PSF_MAP_SOURCE_CONTAINER_H
#define __PSF_MAP_SOURCE_CONTAINER_H

#include "../Core/SharedLibraryExportMacros.h"
#include "MapSource.h"
#include "../IO/H5IODataTree.h"
#include "../IO/OutputArray.h"
#include <vector>

namespace PSF {

    ///A container full of MapSource objects.
    class LIB_PUBLIC MapSourceContainer : public std::vector<MapSource> {
    public:
        ///Create an empty container.
        MapSourceContainer() {}

        ///Initialize the container from an H5IODataTree.
        MapSourceContainer(
            ///The IO tree containing the sources, the configuration for the PSF
            ///map etc.
            const IO::H5IODataTree &data_tree,

            ///The number of apertures to set for the new MapSource objects.
            unsigned num_apertures
        );

        ///\brief All quantities needed to construct the source list from an 
        ///I/O data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    };

} //End PSF namespace.

#endif
