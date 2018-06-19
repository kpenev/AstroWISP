/**\file
 *
 * \brief Define the C-interface functions for the FitPSF library.
 *
 * \ingroup FitPSF
 */

#include "CInterface.h"
#include "Config.h"
#include <cstdarg>
#include <string>
#include <sstream>
#include <iostream>

FittingConfiguration *create_psffit_configuration()
{
    FittingConfiguration* result = reinterpret_cast<FittingConfiguration*>(
        new FitPSF::Config(0, NULL)
    );
#ifdef DEBUG
    std::cerr << "craeted config at: " << result << std::endl;
#endif
    return result;
}

void destroy_psffit_configuration(
    FittingConfiguration *configuration
)
{
    std::cerr << "destroying config at: " << configuration << std::endl;
    delete reinterpret_cast<FitPSF::Config*>(configuration);
}

void update_psffit_configuration(FittingConfiguration *target_configuration,
                                 bool prffit,
                                 ...)
{
#ifdef DEBUG
    std::cerr << "updating config at: " << target_configuration << std::endl;
#endif

    FitPSF::Config *configuration =
        reinterpret_cast<FitPSF::Config*>(target_configuration);

    char fit_mode[] = "fitpsf",
         config_file_option[] = "--config-file",
         empty_str[] = "";
    if(prffit) fit_mode[4]='r';
    char **argv_like = new char*[3];
    argv_like[0] = fit_mode;
    argv_like[1] = config_file_option;
    argv_like[2] = empty_str;

#ifdef DEBUG
    std::cerr << "Pretending executable is: " << argv_like[0] << std::endl;
#endif

    configuration->parse(3, argv_like);
    delete[] argv_like;

    opt::options_description config_file_options;

    config_file_options.add(configuration->cmdline_config_options())
                       .add(configuration->hidden_options());

    va_list arg_list;
    va_start(arg_list, prffit);

    for(
        char *param_name = va_arg(arg_list, char*);
        param_name[0] != '\0';
        param_name = va_arg(arg_list, char*)
    ) {
        char *param_value = va_arg(arg_list, char*);

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
            *configuration
        );
    }

    va_end(arg_list);

    opt::notify(*configuration);
}
