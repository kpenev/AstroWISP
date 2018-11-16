#include "CommandLineUtil.h"

namespace PSF {

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  Grid*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        const std::string &grid_string =
            opt::validators::get_single_string(option_strings);
        value = boost::any(IO::parse_grid_string(grid_string));
    }

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  ModelType*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        std::string psf_model_string =
            opt::validators::get_single_string(option_strings);
        std::transform(psf_model_string.begin(),
                       psf_model_string.end(),
                       psf_model_string.begin(),
                       ::tolower);
        if(psf_model_string == "sdk") value = boost::any(SDK);
        else if(psf_model_string == "bicubic") value = boost::any(BICUBIC);
        else if(psf_model_string == "zero") value = boost::any(ZERO);
        else throw opt::validation_error(
            opt::validation_error::invalid_option_value
        );
    }

} //End PSF namespace.
