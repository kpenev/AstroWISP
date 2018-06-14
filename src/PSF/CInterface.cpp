/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup PSF
 */

#include "CInterface.h"

void expand_term_expression(char *expansion_term_expression,
                            char ***term_list,
                            unsigned *num_terms)
{
    using boost::spirit::ascii::space;
    using boost::spirit::qi::phrase_parse;

    typedef std::list< std::string > GeneratorResultType;
    GeneratorResultType terms;
    PSF::TermGenerator::Grammar< std::string::const_iterator,
        GeneratorResultType > generator_grammar;

    std::string string_expression(expansion_term_expression);

    std::string::const_iterator parse_start = string_expression.begin(),
                                parse_end = string_expression.end();

    bool parser_status = phrase_parse(parse_start,
                                      parse_end,
                                      generator_grammar,
                                      space,
                                      terms);

    if ( !parser_status || parse_start != parse_end ) {
        std::ostringstream message;
        message << "PSF terms parsing failed at "
                << std::string(parse_start, parse_end);
        throw Error::CommandLine(message.str());
    }

    *num_terms = terms.size();

    *term_list = new char*[*num_terms];
    unsigned term_index = 0;
    for(
        GeneratorResultType::const_iterator term_iter = terms.begin();
        term_iter != terms.end();
        ++term_iter, ++term_index
    ) {
        (*term_list)[term_index] = new char[term_iter->size() + 1];
        std::strncpy((*term_list)[term_index],
                     term_iter->c_str(),
                     term_iter->size() + 1);
    }
}

void free_term_list(char **term_list, unsigned num_terms)
{
    for(unsigned i = 0; i < num_terms; ++i)
        delete[] term_list[i];

    delete[] term_list;
}

void evaluate_terms(char **term_expressions,
                    unsigned num_terms,
                    char **variable_names,
                    double **variable_values,
                    unsigned num_variables,
                    unsigned num_sources,
                    double *result)
{
    PSF::MapVarListType variables;
    for(unsigned var_index = 0; var_index < num_variables; ++var_index) {
        variables.push_back(
            PSF::MapVariableType(
                variable_names[var_index],
                std::valarray<double>(variable_values[var_index],
                                      num_sources)
            )
        );
    }
    std::vector<PSF::TermValarray> expansion_term_values(num_terms);
    std::vector<std::string> string_expressions(num_terms);
    for(unsigned term_index = 0; term_index < num_terms; ++term_index) 
        string_expressions[term_index] = term_expressions[term_index];
    evaluate_terms(string_expressions.begin(),
                   string_expressions.end(),
                   variables.begin(),
                   variables.end(),
                   expansion_term_values);

    for(unsigned source_index = 0; source_index < num_sources; ++source_index)
        for(unsigned term_index = 0; term_index < num_terms; ++term_index)
                result[source_index * num_terms + term_index] =
                    expansion_term_values[term_index][source_index];
}
