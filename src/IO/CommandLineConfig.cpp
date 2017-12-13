/**\file
 *
 * \brief Define some of the methods of the CommandLineConfig class.
 *
 * \ingroup IO
 */

#include "CommandLineConfig.h"
#include "SubPixHDF5File.h"

namespace IO {

    void CommandLineConfig::parse(int argc, char **argv)
    {
        describe_options();
        _cmdline_config.add_options()
            (
                "io.hdf5_structure",
                opt::value<std::string>(),
                "An xml file defining the structure of the output HDF5 "
                "files."
            )
            (
                "psf.sdk.max-exp-coef",
                opt::value<double>(),
                "Specifies the maximum value of any term appearing in an "
                "exponent when calculating PSF integrals. Larger values "
                "typically result in faster code, but could lead to extreme "
                "numerical round-off errors."
            )
            (
                "psf.sdk.abs-int-precision",
                opt::value<double>(),
                "Absolute precision up to which integrals of SDK PSFs should"
                " be calculated by default (aperture photometry tunes those "
                "for fast processing)."
            )
            (
                "psf.sdk.rel-int-precision",
                opt::value<double>(),
                "Relative precision up to which integrals of SDK PSFs should"
                " be calculated by default (aperture photometry tunes those "
                "for fast processing)."
            );
        _cmdline_only.add_options()
            (
                "config-file,c",
                opt::value<std::string>()->default_value("subpix.cfg"),
                "The file to read configuration from for options not "
                "specified on the command line."
            )
            (
                "help,h",
                "Print this help and exit"
            );
        opt::options_description cmdline_options,
            config_file_options,
            visible_options;
        cmdline_options.add(_cmdline_config)
            .add(_cmdline_only)
            .add(_hidden);
        config_file_options.add(_cmdline_config)
            .add(_hidden);
        visible_options.add(_cmdline_config)
            .add(_cmdline_only);

        opt::store(
            opt::command_line_parser(argc, argv)
            .options(cmdline_options)
            .positional(_positional)
            .run(),
            *this
        );
        opt::notify(*this);
        if(count("help")) {
            std::cout << usage_help(argv[0]) << visible_options << std::endl;
            return;
        }

        std::ifstream config_stream(
            operator[]("config-file").as<std::string>().c_str()
        );
        if(!config_stream) throw Error::CommandLine(
            "Failed to open config file: "
            +
            operator[]("config-file").as<std::string>()
        );
        opt::store(
            opt::parse_config_file(config_stream,
                                   config_file_options,
                                   true),
            *this
        );
        opt::notify(*this);
        SubPixHDF5File::configure(
            operator[]("io.hdf5_structure").as<std::string>().c_str()
        );

        __executable = argv[0];
        if(__executable == "") __executable = "fitpsf";
        __executable.erase(0, __executable.find_last_of('/') + 1);

        __parsed_ok = true;
    }

} //End IO namespace.
