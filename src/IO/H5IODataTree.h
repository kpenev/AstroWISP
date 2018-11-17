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
#include <H5Cpp.h>
#include <boost/property_tree/ptree.hpp>
#include <list>
#include <valarray>
#include <vector>

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

#ifdef VERBOSE_DEBUG
        ~H5IODataTree() {
            std::cerr << "Destroying H5IODataTree at " << this << std::endl;
        }
#endif
    }; //End H5IODataTree class.

    ///\brief List the entries in an IO tree flagging which are filled and which
    ///are empty.
    LIB_PUBLIC std::ostream &operator<<(
        ///The stream to print to.
        std::ostream &os,
        
        ///The tree to report on.
        const IOTreeBase &tree
    );

} //End IO namespace.

#endif
