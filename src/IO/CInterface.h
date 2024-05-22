/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * IO library.
 *
 * \ingroup IO
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "parse_hat_mask.h"

extern "C" {

    ///See Core::MASK_OK
    LIB_PUBLIC extern const char MASK_OK;

    ///See Core::MASK_CLEAR
    LIB_PUBLIC extern const char MASK_CLEAR;


    ///See Core::MASK_FAULT
    LIB_PUBLIC extern const char MASK_FAULT;

    ///See Core::MASK_HOT
    LIB_PUBLIC extern const char MASK_HOT;

    ///See Core::MASK_COSMIC
    LIB_PUBLIC extern const char MASK_COSMIC;

    ///See Core::MASK_OUTER
    LIB_PUBLIC extern const char MASK_OUTER;

    ///See Core::MASK_OVERSATURATED
    LIB_PUBLIC extern const char MASK_OVERSATURATED;

    ///See Core::MASK_LEAKED
    LIB_PUBLIC extern const char MASK_LEAKED;

    ///See Core::MASK_SATURATED
    LIB_PUBLIC extern const char MASK_SATURATED;

    ///See Core::MASK_INTERPOLATED
    LIB_PUBLIC extern const char MASK_INTERPOLATED;

    ///See Core::MASK_BAD
    LIB_PUBLIC extern const char MASK_BAD;

    ///See Core::MASK_ALL
    LIB_PUBLIC extern const char MASK_ALL;

    ///See Core::MASK_NAN
    LIB_PUBLIC extern const char MASK_NAN;

    ///Opaque struct to cast to/from IO::H5IODataTree.
    struct LIB_PUBLIC H5IODataTree;

    ///\brief C-binding alias for IO::parse_hat_mask.
    ///
    ///See IO::parse_hat_mask() for more details.
    LIB_PUBLIC void parse_hat_mask(const char *mask_string,
                                   long x_resolution,
                                   long y_resolution,
                                   char *mask);

    ///Create a tree to hold the results of AstroWISP processing.
    LIB_PUBLIC H5IODataTree *create_result_tree(
        ///The configuration of the tool being used. Should point to an instance
        ///of some sub-class of IO::CommandLineConfig
        void *configuration,

        ///Information about the version of the tool(s) used.
        char *version_info
    );

    ///\brief Free the memory held by a result tree previously created by
    ///create_result_tree()
    LIB_PUBLIC void destroy_result_tree(
        ///The tree to destroy.
        H5IODataTree *tree
    );

    ///\brief Query the result tree for a value (NULL if value is undefined).
    ///
    ///If the value requested is found and not empty, the return value is true,
    ///otherwise, false and result is not touched.
    LIB_PUBLIC bool query_result_tree(
        ///The tree to extract the quantity from.
        H5IODataTree *tree,

        ///What quantity to get from the tree
        const char *quantity,

        ///Data type of quantity to expect.
        ///Examples:
        ///    "int" for a single integer value
        ///    "[int]" for an array of integers
        ///    "double:double" for a pair of doubles
        ///    "str" for a string entry
        const char *format,

        ///Destination to fill with the value of the quantity. For all formats
        ///except "str", result must be pointer to the correct format and have
        ///sufficient space pre-allocated. For strings, result should be char**,
        ///and new memory is allocated for the string, so the caller is
        ///responsible for free-ing that.
        void *result
    );

    ///\brief Add or overwrite a quantity in/to the result tree.
    ///
    ///If a new quantity was created, the return value is true, if a quantity
    ///was overwritten it is false.
    LIB_PUBLIC void update_result_tree(
        ///Which quantity to add/overwrite.
        const char *quantity,

        ///The value(s) to set for the quantity.
        void *value,

        ///Data type of quantity. See same name argument of query_result_tree()
        ///for details.
        const char *format,

        ///The number of entries in the quantity
        unsigned length,

        ///The tree to update.
        H5IODataTree *tree
    );

    ///\brief Get the PSF map variables for a given image index from the given
    ///result tree.
    ///
    ///If the correspnoding entry in the tree is empty or non-existent, return
    ///false and do not touch the column_data argument.
    LIB_PUBLIC bool get_psf_map_variables(
        ///The tree to extract the varibales from.
        H5IODataTree *output_data_tree,

        ///The image index for which to extract the variables.
        unsigned image_index,

        ///The location to fill with the values of the variables. The values of
        ///each varibale are consecutive in memory. All required storage must
        ///already be allocated.
        double *column_data
    );

    LIB_PUBLIC void set_psf_map_variables(
        ///The names of the variables required by the map.
        char **map_variable_names,

        ///The values of the variables required by the map for each source.
        double *map_variable_values,

        ///The number of map variables
        unsigned num_map_variables,

        ///The number of sources.
        unsigned num_sources,

        ///The index of the image to set the PSF map for.
        unsigned image_index,

        ///The tree to update.
        H5IODataTree *tree
    );

    ///\brief List the names of all quantities that currently have a value in a
    ///data tree.
    ///
    ///The return value is the number of defined quantities found.
    LIB_PUBLIC unsigned list_tree_quantities(
        ///The tree to extract the quantity from.
        H5IODataTree *tree,

        ///A newly allocated array containing all the quantity names. The caller
        ///is responsible for free-ing it when no longer needed.
        char ***quantities
    );

    LIB_PUBLIC void export_free(void *);

} //End Extern "C".
