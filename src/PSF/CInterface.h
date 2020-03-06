/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * PSF library.
 *
 * \ingroup PSF
 */

#include "TermCalculator.h"
#include "Typedefs.h"
#include "../IO/CInterface.h"
#include <string>

extern "C" {
    ///Opaque struct to cast to/from PSF::PiecewiseBicubicMap.
    struct LIB_PUBLIC PiecewiseBicubicPSFMap;

    ///Opaque struct to cast to/from PSF::PiecewiseBicubic.
    struct LIB_PUBLIC PiecewiseBicubicPSF;

    ///\brief Expand the given expression into a list of expressions, each
    ///defining a single PSF expansion term.
    LIB_PUBLIC void expand_term_expression(
        ///The expression to parse. See TermGenerator::Grammar for a description
        ///of the syntax.
        char *expansion_term_expression,

        ///The individual terms in the expression. Allocates new memory.
        ///Caller must release allocated memory when no longer required.
        char ***term_list,

        ///On exit, set to the number of terms found in the term expression.
        unsigned *num_terms
    );

    ///Free the memory allocated by a previous call to expand_term_expression.
    LIB_PUBLIC void free_term_list(
        ///The list of terms created by expand_term_expression() to free.
        char **term_list,

        ///The number of terms in  term_list.
        unsigned num_terms
    );

    ///\brief Evaluate the given terms for the given list of sources.
    LIB_PUBLIC void evaluate_terms(
        ///The list of terms to evaluate. Usually created by
        ///expand_term_expression().
        char **term_list,

        ///The number of terms.
        unsigned num_terms,

        ///The names of the variables used by the expressions in terms.
        char **variable_names,

        ///The values of the variables used by the expressions in terms for each
        ///source. The values of the first variable for all sources should go
        //first, followed by the value of the second variable for all sources
        ///etc.
        double **variable_values,

        ///The number of different variables used.
        unsigned num_variables,

        ///The number of sources to evaluate the terms of.
        unsigned num_sources,

        ///The location to fill with the evaluated terms. Must already be
        ///allocated with a size of num_variables * num_sources. The values of
        ///all terms for the first source are at the beginning of the array,
        ///followed by the value of all terms for the second source etc.
        double *result
    );

    ///Create an instance of PSF::PiecewiseBicubicMap from a recent fit.
    LIB_PUBLIC PiecewiseBicubicPSFMap *create_piecewise_bicubic_psf_map(
        ///The result tree returned by a PSF/PRF fit.
        H5IODataTree *fit_result_tree
    );

    ///\brief Free the memory held by a PSF map previously created by
    ///create_piecewise_bicubic_psf_map()
    LIB_PUBLIC void destroy_piecewise_bicubic_psf_map(
        ///The PSF map to destroy.
        PiecewiseBicubicPSFMap *map
    );

    ///\brief Return a newly allocated PSF/PRF per the given map at the given
    ///location.
    LIB_PUBLIC PiecewiseBicubicPSF *evaluate_piecewise_bicubic_psf_map(
        ///The PSF/PRF map to evaluate.
        PiecewiseBicubicPSFMap *map,

        ///The values of the terms at which to evaluate the map. Usually created
        ///by evaluate_terms for a single source.
        double *term_values
    );

    ///\brief De-allocate a PSF/PRF allocated using
    ///evaluate_piecewise_bicubic_psf_map()
    LIB_PUBLIC void destroy_piecewise_bicubic_psf(
        ///The PSF to delete.
        PiecewiseBicubicPSF *psf
    );

    ///Evaluate a PSF/PRF at a collection of offsets from the source center.
    LIB_PUBLIC void evaluate_piecewise_bicubic_psf(
        ///The PSF/PRF to evaluate.
        PiecewiseBicubicPSF *psf,

        ///The offsets from the source center in the x direction of the points to
        ///evaluate the PSF at. Must have a size equal to num_points.
        double *x_offsets,

        ///The offsets from the source center in the y direction of the points to
        ///evaluate the PSF at. Must have a size equal to num_points.
        double *y_offsets,

        ///The number of locations we are evaluating the PSF/PRF at.
        unsigned num_points,

        ///The location to fill with the values of the PSF/PRF. Must already be
        ///allocated with a size of num_points.
        double *result
    );

    ///\brief See PSF::PSF::integrate, but calculates multiple integrals.
    LIB_PUBLIC void integrate_piecewise_bicubic_psf(
        ///The PSF to integrate.
        PiecewiseBicubicPSF *psf,

        ///The x coordinates of the centers of the rectangles to integrate
        ///over.
        double *center_x,

        ///The y coordinates of the centers of the rectangles to integrate
        ///over.
        double *center_y,

        ///The widths of the rectangles.
        double *dx,

        ///The heights of the rectangles.
        double *dy,

        ///The radii of the circles. For zero entries, the integral is over the
        ///full rectangle.
        double *circle_radii,

        ///The number of integrations requested.
        unsigned num_integrals,

        ///The location to fill with the values of the calculated integrals.
        ///Must already be allocated with a size of num_integrals.
        double *result
    );

};//End extern "C".
