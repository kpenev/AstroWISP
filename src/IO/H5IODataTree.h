/**\file
 *
 * \brief A structure to hold any output data produced by subpixtools.
 *
 * \ingroup IO
 */

#ifndef __H5_OUTPUT_DATA_TREE
#define __H5_OUTPUT_DATA_TREE

#include "../Core/SharedLibraryExportMacros.h"
#include "CommandLineConfig.h"
#include "TranslateToAny.h"
#include "parse_grid.h"
#include "../Background/Annulus.h"
#include "../Core/PhotColumns.h"
#include "../Core/Error.h"
#include "../PSF/Grid.h"
#include "../PSF/Typedefs.h"
#include "../Core/Typedefs.h"
#include "Eigen/Dense"
#include <boost/property_tree/ptree.hpp>
#include <list>
#include <valarray>
#include <vector>
#include <cmath>

// namespace boost {
//   template<class T,
//     typename std::enable_if<std::is_same<T, any>{}, bool>::type =true
//   >
//   bool operator == (const T& lhs, const T& rhs){
//     throw Error::NotImplemented("Comparison of boost any is undefined");
//   }
// }

namespace boost {
  bool operator==(const any& lhs, const any& rhs);
}

namespace IO {

    namespace opt = boost::program_options;


    ///Convenience alias for the boost base class for the IO tree.
    typedef boost::property_tree::basic_ptree<std::string, boost::any>
        IOTreeBase;

    ///\brief A property tree to hold all datasets, attributes and links to
    ///output.
    class LIB_PUBLIC H5IODataTree : public IOTreeBase {
    private:
        ///Tags for the various tools that can fill this tree with data.
        enum TOOL {
            FITPSF,
            FITPRF,
            SUBPIXPHOT,
            FITSUBPIX
        };

        ///The that was used to generate the data in the tree.
        TOOL __tool;

        ///The first part of the key for elements corresponding to the given tool
        std::string __prefix;

        ///The PSF model used (for PSF fitting only).
        std::string __psf_model;

        ///A list of allocated vectors containing string to destroy when the
        ///tree is destroyed;
        std::list< std::vector<std::string> *> __strings_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<int>* > __ints_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<unsigned>* > __uints_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<long>* > __longs_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<unsigned long>* > __ulongs_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<short>* > __shorts_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<unsigned short>* > __ushorts_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<char>* > __chars_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<unsigned char>* > __uchars_to_destroy;

        ///See __strings_to_destroy;
        std::list< std::vector<bool>* > __bools_to_destroy;

        ///See __strings_to_destroy;
        std::list< Eigen::VectorXd* > __doubles_to_destroy;

        ///Destroy one of the __*_to_destroy lists.
        template<typename ARRAY_TYPE>
            void destroy_allocated(std::list< ARRAY_TYPE* > &target);

        ///Add values from a C-style array to the tree.
        template<typename UNIT_TYPE>
            void add_1d_entry(
                ///C-style array of the values to add.
                UNIT_TYPE *value,

                ///The number of entries value
                unsigned length,

                ///The path within the tree to add/update
                const std::string &quantity,

                ///A list to add the newly allocated data to. Should be one of
                //the __*_to_destroy lists.
                std::list< std::vector<UNIT_TYPE>* > &destroy_list
            );

        ///\brief Add values from a C-style double array to the tree.
        void add_1d_entry(
            ///C-style array of the values to add.
            double *value,

            ///The number of entries value
            unsigned length,

            ///The path within the tree to add/update
            const std::string &quantity
        );

        ///Add values from a C-style array of strings to the tree.
        void add_1d_string_entry(
            ///The C-style array of stings to add
            char **value,

            ///How many stings are in the array
            unsigned length,

            ///the path within the tree to add/update.
            const std::string &path
        );

        ///\brief Prepares the tree for the specific tool used.
        void initialize_command_line(
            ///The number of command line tokens.
            int argc,

            ///The command line tokens.
            char** argv,

            ///The executable invoked (no path).
            const std::string &executable,

            ///The version of the tool used.
            const std::string &version
        );

        ///Decides what to do with a single options entry for psf fitting.
        void process_psffit_option(
            ///The key this option is identified by.
            const std::string &key,

            ///The value of the option.
            const opt::variable_value &value
        );

        ///Decides what to do with a single options entry for psf fitting.
        void process_subpixphot_option(
            ///The key this option is identified by.
            const std::string &key,

            ///The value of the option.
            const opt::variable_value &value
        );

    public:
        ///Creates an empty tree.
        H5IODataTree() {}

        ///Fills all command line information.
        H5IODataTree(
            ///The number of command line tokens.
            int argc,

            ///The command line tokens.
            char** argv,

            ///The version of the tool used.
            const std::string &version,

            ///The parsed command line options.
            const CommandLineConfig& options
        )
        {
            initialize_command_line(argc, argv, options.executable(), version);
            fill_configuration(options);
        }

        ///Fills all attributes defining the configuration from the command line.
        void fill_configuration(
            ///The parsed command line options.
            const boost::program_options::variables_map& options
        );

        ///Add a C-style array of values to the tree.
        void add_c_array(
            ///The path within  the tree to add/update.
            const std::string &quantity,

            ///The beginning of the memory where the values are to be found.
            void *value,

            ///Identifier for the type of values being added. For example 'str'
            ///or 'int'.
            const std::string &format,

            ///How many values are in the array.
            unsigned length
        );

        ~H5IODataTree();
    }; //End H5IODataTree class.

    ///\brief List the entries in an IO tree flagging which are filled and which
    ///are empty.
    LIB_PUBLIC std::ostream &operator<<(
        ///The stream to print to.
        std::ostream &os,

        ///The tree to report on.
        const IOTreeBase &tree
    );

    template<typename ARRAY_TYPE>
        void H5IODataTree::destroy_allocated(std::list< ARRAY_TYPE* > &target)
        {
            for(
                typename std::list< ARRAY_TYPE* >::iterator
                target_i = target.begin();
                target_i != target.end();
                ++target_i
            )
                delete *target_i;
        }

    template<typename UNIT_TYPE>
        void H5IODataTree::add_1d_entry(
            UNIT_TYPE *value,
            unsigned length,
            const std::string &quantity,
            std::list< std::vector<UNIT_TYPE>* > &destroy_list
        )
        {
            std::vector<UNIT_TYPE> *entry=new std::vector<UNIT_TYPE>(
                value,
                value + length
            );
            if(length > 1) {
                put(quantity,
                    entry,
                    IO::TranslateToAny< std::vector<UNIT_TYPE> >());
                destroy_list.push_back(entry);
            } else if(length == 1) {
                put(quantity,
                    (*entry)[0],
                    IO::TranslateToAny< UNIT_TYPE >());
                delete entry;
            } else
                throw Error::InvalidArgument(
                    "add_1d_tree_entry",
                    "Attempting to add zero length dataset to I/O tree!"
                );
        }

} //End IO namespace.

#endif
