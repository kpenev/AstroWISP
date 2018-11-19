/**\file
 *
 * \brief Declarations of the Map class and some related functions.
 */

#ifndef __PSF_MAP_H
#define __PSF_MAP_H

#include "../Core/SharedLibraryExportMacros.h"
#include "PSF.h"
#include "MapSource.h"
#include "../Core/Point.h"
#include "../Core/Error.h"
#include "Eigen/Dense"
#include "../IO/H5IODataTree.h"
#include <set>

namespace PSF {

    ///\brief Interface for working with smooth dependenc of PSF/PRF on source
    ///parameters.
    class LIB_PUBLIC Map {
    private:
        ///See num_terms()
        unsigned __num_terms;
    public:
        ///Create a map with the given number of terms.
        Map(
            ///See num_terms()
            unsigned num_terms=0
        ) : __num_terms(num_terms) {}

        virtual ~Map() {}

        ///The number of terms that the map depends on.
        virtual unsigned num_terms() const {return __num_terms;}

        ///Set the number of terms that the map depends on.
        virtual void set_num_terms(
            ///The number of terms to set.
            unsigned nterms
        )
        {__num_terms = nterms;}

        ///A reference to a dynamically allocated PSF.
        virtual PSF *operator()(
            ///The values of the terms on which the PSF map depends.
            const Eigen::VectorXd &terms,

            ///Background to add to the PSF (assumed constant).
            double background = 0
        ) const =0;

        ///A reference to a dynamically allocated PSF.
        PSF *operator()(
            ///The source whose PSF we want.
            const MapSource &source,

            ///Background to add to the PSF (assumed constant) normalized
            ///in the same way as the backgroundless PSFs produced by
            ///the map.
            double background = 0
        ) const
        {return (*this)(source.expansion_terms(), background);}

        ///\brief All quantities needed to construct the PSF map from an I/O
        ///data tree.
        static const std::set<std::string> &required_data_tree_quantities();
    };

    ///\brief Combine two collections of required data tree quantities into a set.
    ///
    ///Used by inheriting PSF maps to add their required quantities.
    std::set<std::string> combine_required_tree_quantities(
        ///The first iterator from the first set to include.
        std::set<std::string>::const_iterator v1_start,

        ///One past the last iterator from the first set to include.
        std::set<std::string>::const_iterator v1_end,

        ///The first iterator from the second set to include.
        const std::string* v2_start,

        ///One past the last iterator from the second set to include.
        const std::string* v2_end
    );

    ///A simpler version of combine_required_tree_quantities().
    std::set<std::string> combine_required_tree_quantities(
        const std::string* extra_start,
        const std::string* extra_end
    );

} //End PSF namespace.

#endif
