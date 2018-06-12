/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup PSF
 */

#include "CInterface.h"

void expand_term_expression(char *expansion_term_expression,
                            char **term_list,
                            unsigned *num_terms)
{
    using boost::spirit::ascii::space;
    using boost::spirit::qi::phrase_parse;

    typedef std::list< std::string > GeneratorResultType;
    GeneratorResultType terms;
    TermGenerator::Grammar< std::string::const_iterator,
        GeneratorResultType > generator_grammar;

    bool parser_status = phrase_parse(
        expansion_term_expression
        (
            expansion_term_expression
            +
            std::char_traits<char>::length(expansion_term_expression)
        ),
        generator_grammar,
        space,
        terms
    );

    num_terms = terms.size();

    term_list = new char*[num_terms];
    unsigned term_index = 0;
    for(
        GeneratorResultType::const_iterator term_iter = terms.begin();
        term_iter != terms.end();
        ++term_iter, ++term_index
    ) {
        term_list[term_index] = new char[term_iter->size() + 1];
        std::strcpy(term_list[term_index],
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
                    double **variables,
                    unsigned num_variables,
                    unsigned num_sources
                    double *result)
{
}
