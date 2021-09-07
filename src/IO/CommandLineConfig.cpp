/**\file
 *
 * \brief Define some of the methods of the CommandLineConfig class.
 *
 * \ingroup IO
 */

#include "CommandLineConfig.h"
#include <cstdarg>

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

        if(!(operator[]("config-file").as<std::string>().empty())) {
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
        }

        __executable = argv[0];
        if(__executable == "") __executable = "fitpsf";
        __executable.erase(0, __executable.find_last_of('/') + 1);

        __parsed_ok = true;
    }

    void update_configuration(IO::CommandLineConfig &configuration,
                              const std::string &mock_executable,
                              va_list options)
    {
#ifdef DEBUG
        std::cerr << "updating config at: " << &configuration << std::endl;
#endif

        char config_file_option[] = "--config-file",
             empty_str[] = "";
        char **argv_like = new char*[3];
        argv_like[0] = const_cast<char*>(mock_executable.c_str());
        argv_like[1] = config_file_option;
        argv_like[2] = empty_str;

#ifdef DEBUG
        std::cerr << "Pretending executable is: " << argv_like[0] << std::endl;
#endif

        configuration.parse(3, argv_like);
        delete[] argv_like;

        opt::options_description config_file_options;

        config_file_options.add(configuration.cmdline_config_options())
            .add(configuration.hidden_options());


        for(
            char *param_name = va_arg(options, char*);
            param_name[0] != '\0';
            param_name = va_arg(options, char*)
        ) {
            char *param_value = va_arg(options, char*);

            std::stringstream config_stream;
            config_stream << param_name << " = " << param_value << std::endl;
#ifdef DEBUG
            config_stream.seekg(0, std::ios::beg);
            std::cerr << "Setting " << config_stream.str() << std::endl;
#endif
            config_stream.seekg(0, std::ios::beg);
            opt::store(
                opt::parse_config_file(config_stream,
                                       config_file_options,
                                       true),
                configuration
            );
        }

        opt::notify(configuration);
    }

} //End IO namespace.
