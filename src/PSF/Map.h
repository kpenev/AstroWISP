/**\brief
 *
 * \file Declarations of the Map class and some related functions.
 */

#ifndef __PSF_MAP_H
#define __PSF_MAP_H

#include "PSF.h"
#include "MapSource.h"
#include "../Core/Point.h"
#include "../Core/Error.h"
#include "Eigen/Dense"
#include "../IO/H5IODataTree.h"
#include <set>

namespace PSF {

    class Map {
    private:
        ///The number of terms that the map depends on.
        unsigned __num_terms;
    public:
        Map(unsigned num_terms=0) : __num_terms(num_terms) {}

        virtual ~Map() {}

        ///The number of terms that the map depends on.
        virtual unsigned num_terms() const {return __num_terms;}

        ///Set the number of terms that the map depends on.
        virtual void set_num_terms(unsigned nterms)
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

    ///To be used by inheriting PSF maps to add their required quantities.
    std::set<std::string> combine_required_tree_quantities(
        std::set<std::string>::const_iterator v1_start,
        std::set<std::string>::const_iterator v1_end,
        const std::string* v2_start,
        const std::string* v2_end
    );

    ///A simpler version of combine_required_tree_quantities().
    std::set<std::string> combine_required_tree_quantities(
        const std::string* extra_start,
        const std::string* extra_end
    );

} //End PSF namespace.

#endif
