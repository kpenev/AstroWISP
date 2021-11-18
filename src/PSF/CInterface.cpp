/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h
 *
 * \ingroup PSF
 */

#include "CInterface.h"
#include "PiecewiseBicubicMap.h"
#include "PiecewiseBicubic.h"

PiecewiseBicubicPSFMap *create_piecewise_bicubic_psf_map(
    H5IODataTree *fit_result_tree
)
{
    return reinterpret_cast<PiecewiseBicubicPSFMap*>(
        new PSF::PiecewiseBicubicMap(
            *reinterpret_cast<IO::H5IODataTree*>(fit_result_tree)
        )
    );
}

void destroy_piecewise_bicubic_psf_map(PiecewiseBicubicPSFMap *map)
{
    delete reinterpret_cast<PSF::PiecewiseBicubicMap*>(map);
}

PiecewiseBicubicPSF *evaluate_piecewise_bicubic_psf_map(
    PiecewiseBicubicPSFMap *map,
    double *term_values
)
{
    PSF::PiecewiseBicubicMap *real_map =
        reinterpret_cast<PSF::PiecewiseBicubicMap*>(map);
    return reinterpret_cast<PiecewiseBicubicPSF*>(
        real_map->get_psf(
            Eigen::Map<const Eigen::VectorXd>(term_values, real_map->num_terms())
        )
    );
}

void destroy_piecewise_bicubic_psf(PiecewiseBicubicPSF *psf)
{
    delete reinterpret_cast<PSF::PiecewiseBicubic*>(psf);
}

void evaluate_piecewise_bicubic_psf(PiecewiseBicubicPSF *psf,
                                    double *x_offsets,
                                    double *y_offsets,
                                    unsigned num_points,
                                    double *result)
{
    PSF::PiecewiseBicubic *real_psf = reinterpret_cast<PSF::PiecewiseBicubic*>(
        psf
    );
    for(unsigned i = 0; i < num_points; ++i) {
        result[i] = (*real_psf)(x_offsets[i], y_offsets[i]);
    }
}

void integrate_piecewise_bicubic_psf(PiecewiseBicubicPSF *psf,
                                     double *center_x,
                                     double *center_y,
                                     double *dx,
                                     double *dy,
                                     double *circle_radii,
                                     unsigned num_integrals,
                                     double *result)
{
    PSF::PiecewiseBicubic *real_psf = reinterpret_cast<PSF::PiecewiseBicubic*>(
        psf
    );
    for(unsigned i = 0; i < num_integrals; ++i) {
        result[i] = (*real_psf).integrate(center_x[i],
                                          center_y[i],
                                          dx[i],
                                          dy[i],
                                          circle_radii[i]);
    }
}
