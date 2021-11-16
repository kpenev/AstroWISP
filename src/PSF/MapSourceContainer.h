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
    class LIB_LOCAL MapSourceContainer : public std::vector<MapSource> {
    public:
        ///Create an empty container.
        MapSourceContainer() {}

        ///Initialize the container from an H5IODataTree.
        MapSourceContainer(
            ///The IO tree containing the sources, the configuration for the PSF
            ///map etc.
            const IO::H5IODataTree &data_tree,

            ///The number of apertures to set for the new MapSource objects.
            unsigned num_apertures,

            ///The string to add to node names when querying data_tree to select
            ///the entries corresponding to the image being processed within the
            ///output tree. If no split by image is present, use an empty
            ///string.
            const std::string &data_tree_image_id=""
        );

        ///\brief All quantities needed to construct the source list from an
        ///I/O data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    };

} //End PSF namespace.

#endif
