/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * PSF library.
 *
 * \ingroup PSF
 */

#include "TermCalculator.h"
#include <string>

extern "C" {
    ///\brief Expand the given expression into a list of expressions, each
    ///defining a single PSF expansion term.
    void expand_term_expression(
        ///The expression to parse. See TermGenerator::Grammar for a description
        ///of the syntax.
        char *expansion_term_expression,

        ///The individual terms in the expression. Allocates new memory.
        ///Caller must release allocated memory when no longer required.
        char **term_list,

        ///On exit, set to the number of terms found in the term expression.
        unsigned *num_terms
    );

    ///Free the memory allocated by a previous call to expand_term_expression.
    void free_term_list(
        ///The list of terms created by expand_term_expression() to free.
        char **term_list,

        ///The number of terms in  term_list.
        unsigned num_terms
    );

    ///\brief Evaluate the given terms for the given list of sources.
    void evaluate_terms(
        ///The list of terms to evaluate. Usually created by
        ///expand_term_expression().
        char **term_list,

        ///The number of terms.
        unsigned num_terms,

        ///The names of the variables used by the expressions in terms.
        char **variable_names,

        ///The values of the variables used by the expressions in terms for each
        ///source.
        double **variables,

        ///The number of different variables used.
        unsigned num_variables,

        ///The number of sources to evaluate the terms of.
        unsigned num_sources

        ///The location to fill with the evaluated terms. Must already be
        ///allocated with a size of num_variables * num_sources.
        double *result
    );

};//End extern "C".
