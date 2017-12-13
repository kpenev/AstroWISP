/**\file
 * 
 * \brief Declare a class for managing the configuration with which tools are 
 * run.
 *
 * \ingroup IO
 */

#ifndef __COMMAND_LINE_CONFIG_H
#define __COMMAND_LINE_CONFIG_H

#include "CommandLineUtil.h"
#include <boost/program_options.hpp>
#include "../PSF/EllipticalGaussian.h"


namespace opt = boost::program_options;

namespace IO {

    ///A base class for configuration of tools using the command line.
    class CommandLineConfig : public opt::variables_map {
    private:
        ///Whether the parsing was successful.
        bool __parsed_ok;

        ///The name of the executable invoked (without any path information).
        std::string __executable;

        ///Describes the available command line options.
        virtual void describe_options() =0;

        ///\brief The part of the help describing the usage and purpose (no 
        ///options).
        virtual std::string usage_help(
            ///The name of the executable being run.
            const std::string &
        ) const
        {return "";};

    protected:
        opt::options_description
            ///\brief The descriptions of options not included in the help 
            ///message.
            ///
            ///In this case the positional options.
            _hidden,

            ///Options accessible only on the command line.
            _cmdline_only,

            ///Options accessible from both command line and config file.
            _cmdline_config;

        ///Description of the positional command line arguments.
        opt::positional_options_description _positional;

    public:
        CommandLineConfig(
            int argc = 0,
            char **argv = NULL
        )
            : __parsed_ok(false)
        {if(argc > 0) parse(argc, argv);}

        ///Parse the command line.
        void parse(
            ///The number of arguments on the command line
            ///(+1 for the executable)
            int argc,

            ///A C style array of the actual command line arguments.
            char **argv
        );

        ///Did parsing succeed?
        bool proceed() const {return __parsed_ok;}

        ///The executable whose command line is being processed (no path).
        const std::string &executable() const {return __executable;}
    }; //End CommandLineConfig class.

} //End IO namespace.

#endif
