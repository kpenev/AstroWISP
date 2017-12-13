#include "CommandLineUtil.h"

namespace Background {

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  Annulus*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        const std::string &bg_annulus_string = 
            opt::validators::get_single_string(option_strings);
        std::list<double> bgan_values = Core::parse_real_list(
            bg_annulus_string,
            "--bg-annulus",
            2,
            2
        );
        value = boost::any(
            Annulus(bgan_values.front(), bgan_values.back())
        );
    }

} //End Background namespace.
