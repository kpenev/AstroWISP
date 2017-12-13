/**\file
 *
 * \brief Declares a container for holding MapSource objects.
 *
 * \ingroup PSF
 */

#include "MapSource.h"
#include "../IO/H5IODataTree.h"
#include "../IO/OutputArray.h"
#include <vector>

namespace PSF {

    ///A container full of MapSource objects.
    class MapSourceContainer : public std::vector<MapSource> {
    public:
        ///Create an empty container.
        MapSourceContainer() {}

        ///Initialize the container from an H5IODataTree.
        MapSourceContainer(const IO::H5IODataTree &data_tree,
                           unsigned num_apertures);

        ///\brief All quantities needed to construct the source list from an 
        ///I/O data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    };

} //End PSF namespace.
