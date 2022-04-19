/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * PSF library.
 *
 * \ingroup PSF
 */

#include "Typedefs.h"
#include "../IO/CInterface.h"
#include <string>

extern "C" {
    ///Opaque struct to cast to/from PSF::PiecewiseBicubicMap.
    struct LIB_PUBLIC PiecewiseBicubicPSFMap;

    ///Opaque struct to cast to/from PSF::PiecewiseBicubic.
    struct LIB_PUBLIC PiecewiseBicubicPSF;

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
