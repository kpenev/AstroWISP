#include "TermCalculator.h"


namespace PSF {

    unsigned TermValarray::__max_id = 0;

    std::ostream &operator<<(std::ostream &os,
                             const std::valarray<double> &array)
    {
        os << "{";
        for ( size_t i=0; i < array.size(); ++i)
            os << array[i] << ( i < array.size() - 1 ? ", " : "");
        os << "}";
        return os;
    }

#ifdef STANDALONE
    void define_some_variables(MapVarListType &variables,
                               unsigned length = 3)
    {
        variables.push_back(
                MapVariableType("pi", TermValarray(M_PI, length))
        );
        variables.push_back(
                MapVariableType("seq", TermValarray(length))
        );
        for(unsigned i=0; i<length; ++i) variables.back().second[i] = i;
    }

    void exercise_term_calculator()
    {
        MapVarListType variables;
        define_some_variables(variables);
        std::vector<TermValarray> result;
        std::string to_parse;
        do {
            std::cout << "Expression> ";
            std::cout.flush();
            std::getline(std::cin, to_parse);
            if ( to_parse != "" )
                try {
                    evaluate_terms(&to_parse, &to_parse + 1,
                                   variables.begin(), variables.end(),
                                   result);
                    std::cout << result.back() << std::endl;
                } catch(Error::ParsingError &error) {
                    std::cout << error.what() << ":"
                              << error.get_message()
                              << std::endl;
                }
        } while ( to_parse != "" );
    }

    void exercise_term_generator()
    {
        std::string to_parse;
        do {
            std::cout << "Expression> ";
            std::cout.flush();
            std::getline(std::cin, to_parse);
            if ( to_parse != "" ) {
                using boost::spirit::ascii::space;
                using boost::spirit::qi::phrase_parse;
                std::string::const_iterator parse_start=to_parse.begin(),
                                            parse_end=to_parse.end();

                typedef std::list< std::string > ResultType;
                ResultType result;

                TermGenerator::Grammar< std::string::const_iterator,
                                        ResultType > generator_grammar;

                bool parser_status = phrase_parse(parse_start,
                                                  parse_end,
                                                  generator_grammar,
                                                  space,
                                                  result);

                if ( parser_status && parse_start == parse_end ) {
                    std::cout << "Terms:" << std::endl;
                    for (
                            ResultType::const_iterator i=result.begin();
                            i!=result.end();
                            ++i
                    ) std::cout << "    " << *i << std::endl;
                } else
                    std::cout << "Parsing failed at "
                              << std::string(parse_start, parse_end)
                              << std::endl;
            }
        } while ( to_parse != "" );
    }

    void exercise()
    {
        std::string to_parse;
        MapVarListType variables;
        define_some_variables(variables);

        do {
            std::cout << "Expression> ";
            std::cout.flush();
            std::getline(std::cin, to_parse);

            std::vector<TermValarray> result(terms.size());

            if ( to_parse != "" )
                try {
                    evaluate_term_expression(to_parse);
                    for(
                        unsigned term_i = 0;
                        term_i < result.size();
                        ++term_i
                    ) {
                        std::cout << "(" << term_i << "): "
                                  << result[term_i] << std::endl;
                    }
                } catch(Error::ParsingError &error) {
                    std::cout << error.what() << ":"
                              << error.get_message()
                              << std::endl;
                } catch(Error::CommandLine &error) {
                        std::cout << error.get_message();
                }

        } while ( to_parse != "" );
    }

    int main()
    {
        std::cout << "Term generator grammar: \n\n"
                  << TermGenerator::ebnf_definition
                  << std::endl;
        exercise();
    }
#endif

}
