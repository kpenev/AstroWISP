#include "CommandLineUtil.h"

namespace Core {

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  RealList*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        const std::string &real_list_string =
            opt::validators::get_single_string(option_strings);
        try {
            value=boost::any(
                    parse_real_list(real_list_string,
                                    "",
                                    0,
                                    std::numeric_limits<unsigned>::max())
            );
        } catch(const Error::CommandLine &) {
            throw opt::validation_error(
                    opt::validation_error::invalid_option_value
            );
        }
    }

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  StringList*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        const std::string &string_list_string =
            opt::validators::get_single_string(option_strings);
        try {
            value=boost::any(
                    parse_string_list(string_list_string, "", 0,
                                      std::numeric_limits<unsigned>::max())
            );
        } catch(const Error::CommandLine &) {
            throw opt::validation_error(
                opt::validation_error::invalid_option_value
            );
        }
    }

    void validate(boost::any& value,
                  const std::vector<std::string>& option_strings,
                  ColumnList*,
                  int)
    {
        opt::validators::check_first_occurrence(value);
        const std::string &column_list_string =
            opt::validators::get_single_string(option_strings);
        try {
            value = boost::any(parse_column_list(column_list_string,
                                                 1,
                                                 "",
                                                 true));
        } catch(const Error::CommandLine &) {
            throw opt::validation_error(
                    opt::validation_error::invalid_option_value
            );
        }
    }

} //End Core namespace.


