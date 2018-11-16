/**\file
 *
 * \brief The declarations of Map methods and related functions.
 */

#include "Map.h"

namespace PSF {

    const std::set<std::string> &Map::required_data_tree_quantities()
    {
        const std::string required_data_tree_quantities[]={
            "psffit.terms"
        };
        static const std::set<std::string> required_quantities(
                required_data_tree_quantities,
                required_data_tree_quantities
                +
                sizeof(required_data_tree_quantities)
                /
                sizeof(required_data_tree_quantities[0])
        );
        return required_quantities;
    }


    std::set<std::string> combine_required_tree_quantities(
            std::set<std::string>::const_iterator v1_start,
            std::set<std::string>::const_iterator v1_end,
            const std::string* v2_start,
            const std::string* v2_end)
    {
        std::set<std::string> result(v1_start, v1_end);
        result.insert(v2_start, v2_end);
        return result;
    }

    std::set<std::string> combine_required_tree_quantities(
            const std::string* extra_start,
            const std::string* extra_end)
    {
        std::set<std::string> result(
            Map::required_data_tree_quantities()
        );
        result.insert(extra_start, extra_end);
        return result;
    }

}
