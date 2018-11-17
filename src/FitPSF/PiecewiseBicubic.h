/**\file
 *
 * \brief The declarations of the functions related to fitting piecewise
 * bicubic PSFs.
 *
 * The following variable names are used for the quentities defined in the
 * [description of the algorithm](@ref PiecewiseBicubicFitPSF_main_page).
 *
 *   - \f$\mathbf{\tilde{M}^i}\f$ - pix_integral_matrix
 *
 *   - The stack of \f$(\mathbf{\kappa}^i)^T\,\mathbf{\kappa}\f$ -
 *     poly_coef_matrix
 *
 *   - \f$\mathbf{\tilde{\Lambda}}\f$ - symmetrized_pix_integral_matrix
 *
 *   - \f$\mathbf{\Lambda}\f$ - matrix_to_invert
 *
 *   - \f$\mathbf{r}\f$ - rhs
 *
 *   - \f$\mathbf{\tilde{r}}\f$ - modified_rhs
 *
 *   - \f$\mathbf{r'}\f$ - flux_scaled_modified_rhs
 *
 * The last four quantities are re-computed after changing the source
 * amplitudes.
 *
 * \ingroup FitPSF
 */

#ifndef __PIECEWISE_BICUBIC_PSF_FITTING_H
#define __PIECEWISE_BICUBIC_PSF_FITTING_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Common.h"
#include "LinearSource.h"
#include "OverlapGroup.h"
#include "PiecewiseBicubicPSFSmoothing.h"
#include "../PSF/PiecewiseBicubic.h"
#include "../PSF/PiecewiseBicubicMap.h"
#include "../Core/SubPixelCorrectedFlux.h"
#include "../Core/SubPixelMap.h"
#include "../Core/Image.h"
#include "../Core/Typedefs.h"
#include "Eigen/Dense"
#include <ctime>
#include <fstream>
#include <set>
#include <cmath>
#include <vector>
#include <valarray>

namespace FitPSF {

    ///Convenience alias for unmutable iterator over pixels.
    typedef std::list< Pixel<LinearSource>* >::const_iterator
        ConstPixelIter;

    ///Convenience alias for mutable iterator over pixels.
    typedef std::list< Pixel<LinearSource>* >::const_iterator
        PixelIter;

    ///Convenience alias for a list of groups of overlapping sources.
    typedef std::list< OverlapGroup<LinearSource, PSF::PiecewiseBicubic> >
        OverlapGroupList;

    ///\brief Generates the vector of bicubic cell parameter sets
    ///corresponding to single PSF parameter=1 while all others=0.
    void LIB_LOCAL fill_parameter_sets(
        ///The vertical grid boundaries defining the PSF cells.
        const std::vector<double> &x_grid,

        ///The horizontal grid boundaries defining the PSF cells.
        const std::vector<double> &y_grid,

        ///The location to fill.
        std::vector< std::valarray<double> > &parameter_sets
    );

    ///\brief Selects basis vectors in which the PSF will be decomposed, as
    ///well as the full PSF integral vector.
    void LIB_LOCAL select_basis_vectors(
        ///The sets of parameters returned by fill_parameter_sets().
        const std::vector< std::valarray<double> > &parameter_sets,

        ///A PSF with the correct grid of cells.
        const PSF::PiecewiseBicubic &psf,

        ///The matrix to fill. The first column will be \f$\mathbf{I}/I^2\f$,
        ///and subsequent columns will be the basis vectors.
        Eigen::MatrixXd &basis
    );

    ///\brief Generates parameter sets which correspond to the basis vectors
    ///(and full PSF integral vector) instead of to single parameters.
    void LIB_LOCAL fill_basis_paramater_sets(
        ///The sets of parameters returned by fill_parameter_sets().
        const std::vector< std::valarray<double> > &parameter_sets,

        ///The matrix of basis vectors filled by select_basis_vectors.
        Eigen::MatrixXd &basis,

        ///The destination to fill with the new parameter sets. Should have
        ///the correct size on input, but the individual entries are
        ///automatically resized.
        std::vector< std::valarray<double> > &basis_parameter_sets
    );

    ///Fills in the 3 output matrices and the output vector as appropriate.
    void LIB_LOCAL prepare_linear_regression(
        ///The input list of sources. Modified because of pixel iterations.
        LinearSourceList &fit_sources,

        ///The sets of parameters returned by fill_parameter_sets().
        const std::vector< std::valarray<double> > &parameter_sets,

        ///The vector of background excess values of the source pixels.
        const Eigen::VectorXd &pixel_excesses,

        ///A stack of matrices, one for each source, with each row consisting
        ///of the integrals over the same source pixel but with only one
        ///basis vector contributing (\f$\mathbf{\tilde{M}^i}\f$ in the
        ///documentation). Should have the correct dimensions on input.
        Eigen::MatrixXd &pixel_integral_matrix,

        /// \f$\mathbf{\Delta}\f$ in the documentation. Should already have
        //the correct size.
        Eigen::VectorXd &rhs_offset,

        ///A stack of \f$(\mathbf{\tilde{M}^i})^T\mathbf{\tilde{M}^i}\f$
        ///(\f$\mathbf{\tilde{\Lambda}}\f$ in the documentation). Should have
        //the correct dimensions on input.
        Eigen::MatrixXd &symmetrized_pix_integral_matrix,

        ///Derived from pixel_excesses by applying
        /// \f$\left(\mathbf{\tilde{M}^i}\right)^T\f$ on each
        ///source's pixels. Should have the correct dimensions on input.
        Eigen::VectorXd &modified_pixel_excesses,

        /// \f$\mathbf{\tilde{\Delta}}\f$ in the documentation.
        Eigen::VectorXd &modified_rhs_offset
    );

    ///\brief Fills a stack of matrices, containing the outer product of the
    ///PSF terms for each source.
    void LIB_LOCAL fill_poly_coef_matrix(
        ///The input list of sources.
        const LinearSourceList &fit_sources,

        ///The location to fill with the result. Should have the correct
        ///dimensions on input.
        Eigen::MatrixXd &poly_coef_matrix
    );

    ///\brief Fills the matrix defining the equations of the linear regression
    ///for the full psf expansion.
    void LIB_LOCAL fill_matrix_to_invert(
        ///The input list of sources.
        const LinearSourceList &fit_sources,

        ///See the argument with the same name in prepare_linear_regression.
        const Eigen::MatrixXd &symmetrized_pix_integral_matrix,

        ///The result produced by fill_poly_coef_matrix.
        const Eigen::MatrixXd &poly_coef_matrix,

        ///The location to fill with the result. Should have the correct
        //dimensions and be all zero on input.
        Eigen::MatrixXd &matrix_to_invert
    );

    ///\brief Creates the RHS vector for the PSF fitting from the pixel
    ///values of the fit sources.
    void LIB_LOCAL fill_pixel_excesses(
        ///The input list of sources. Modified due to pixel iterations.
        LinearSourceList &fit_sources,

        ///The vector to fill. Should already have the correct size.
        Eigen::VectorXd &rhs,

        ///If passed, only the first source in the list is used. This is
        ///useful when defining the RHS for groups of overlapping sources,
        ///which all share a common set of image pixels
        bool single_source = false
    );

    ///\brief Fills a vector with the final RHS that should be used for the
    ///linear regression of the PSF parameter expansion fit.
    void LIB_LOCAL fill_flux_scaled_modified_rhs(
        ///The input list of sources.
        const LinearSourceList &fit_sources,

        ///See the argument with the same name in prepare_linear_regression.
        const Eigen::VectorXd &modified_pixel_excesses,

        ///See the argument with the same name in prepare_linear_regression.
        const Eigen::VectorXd &modified_rhs_offset,

        ///The location to fill with the result. Should have the correct
        ///dimensions and be zero on input.
        Eigen::VectorXd &flux_scaled_modified_rhs
    );

    ///\brief Provides initial guesses for the amplitudes of the input
    ///sources.
    ///
    ///Amplitudes are estimated by summing up all pixel values.
    void LIB_LOCAL estimate_initial_amplitudes(
        ///The input list of sources. On output has the amplitudes set to
        ///their estimated values.
        LinearSourceList &fit_sources,

        ///The gain of the input image.
        double gain
    );

    ///\brief Provides initial guesses for the amplitudes of the input
    ///sources.
    ///
    ///Amplitudes are estimated by aperture photometry with a flat PSF but the
    ///actual sub-pixel map.
    void LIB_LOCAL estimate_initial_amplitudes(
        ///The input list of sources. On output has the amplitudes set to
        ///their estimated values.
        LinearSourceList &fit_sources,

        ///The sub-pixel sensitivity map.
        const Core::SubPixelMap &subpix_map,

        ///The image on which PSF fitting is performed.
        const Core::Image<double> &observed_image,

        ///The gain of the input image.
        double gain,

        ///The aperture to use.
        double aperture
    );

    ///\brief Updates the fluxes of the fit sources with new best estimates
    ///derived by fitting each source's background excesses.
    ///
    ///\return the root sum square change of the fluxes.
    double LIB_LOCAL update_fluxes(
        ///The list of sources for which to update the flux.
        LinearSourceList &fit_sources,

        ///The current best fit polynomial expansion coefficients for the PSF
        ///parameters.
        const Eigen::VectorXd &best_fit
    );

    ///\brief Performs a single PSF parameter fit, amplitude fit iteration.
    ///
    ///\return the sum root square change of the source fluxes.
    double LIB_LOCAL fit_piecewise_bicubic_psf_step(
        ///The list of sources to fit. On output their fluxes are updated to
        ///the new best fit value.
        LinearSourceList &fit_sources,

        ///See prepare_linear_regression.
        const Eigen::MatrixXd &symmetrized_pix_integral_matrix,

        ///See fill_poly_coef_matrix.
        const Eigen::MatrixXd &poly_coef_matrix,

        ///See prepare_linear_regression.
        const Eigen::MatrixXd &pixel_integral_matrix,

        ///See prepare_linear_regression.
        const Eigen::VectorXd &pixel_excesses,

        ///See prepare_linear_regression.
        const Eigen::VectorXd &rhs_offset,

        ///See prepare_linear_regression.
        const Eigen::VectorXd &modified_pixel_excesses,

        ///See prepare_linear_regression.
        const Eigen::VectorXd &modified_rhs_offset,

        ///Correction to the final matrix enforcing smoothing.
        const PiecewiseBicubicPSFSmoothing &smoothing,

        ///Updated to a new estimate of the polynomial expansion coefficients
        ///for the PSF parameters. Should have the correct dimensions on
        ///input.
        Eigen::VectorXd &best_fit
    );

    ///\brief Counts the total number of shape fitting pixels in all the
    ///given sources.
    size_t LIB_LOCAL count_pixels(
        ///The list of sources whose pixel counts should be summed up.
        const LinearSourceList &fit_sources,

        ///If non-NULL the location pointed to by this argument is set to
        ///the number of shape fit pixels in the source with the largest
        ///number of such.
        size_t *max_source_pixels = NULL
    );

    ///\brief Fits for the polynomials giving the values, x, y and xy
    ///derivatives describing a piecewise bicubic PSF smoothly varying over
    ///the image.
    bool LIB_LOCAL fit_piecewise_bicubic_psf(
        ///The input list of sources (need to be modified to perform the
        ///fit) and on output has the amplitudes set correctly for the best
        ///fit.
        LinearSourceList &fit_sources,

        ///The list of sources dropped from the fit (need to be modified to
        ///set the amplitudes correctly for the best fit).
        LinearSourceList &dropped_sources,

        ///The gain of the input image.
        double gain,

        ///The vertical grid boundaries defining the PSF cells.
        const std::vector<double> &x_grid,

        ///The horizontal grid boundaries defining the PSF cells.
        const std::vector<double> &y_grid,

        ///The sub-pixel sensitivity map.
        const Core::SubPixelMap &subpix_map,

        ///The maximum sum squared of amplitude changes allowed after a
        ///PSF/Amplitude fit iteration in order to consider the solution
        ///converged.
        double max_abs_amplitude_change,

        ///The maximum fraction the sum squared of the amplitude changes can
        ///be of the length of the vector of amplitudes after a PSF/Amplitude
        ///fit iteration in order to consider the solution converged.
        double max_rel_amplitude_change,

        ///Sources with reduced \f$\chi^2\f$ bigger than this are removed
        ///from the fit.
        double max_chi2,

        ///Pixels whose residuals after fitting is too large are not included
        ///in the fit. This parameter determines how large the residual
        ///cut-off is. See the help of --pix-rej command line option for
        ///details.
        double pixel_rejection,

        ///The minimum rate of convergence to require before giving up.
        double min_convergence_rate,

        ///The maximum number of iterations before giving up (use a negative
        ///value to disable this limit).
        int max_iterations,

        ///Smoothing penalty (see help for details).
        double smoothing_penalty,

        ///On output contains the best fit polynomial coefficients. Resized
        ///as necessary.
        Eigen::VectorXd &best_fit_poly_coef
    );

    ///\brief Same as above, but uses an initial PSF map to derive initial
    ///flux estimates.
    bool LIB_LOCAL fit_piecewise_bicubic_psf(
        ///The input list of sources (need to be modified to perform the
        ///fit) and on output has the amplitudes set correctly for the best
        ///fit.
        LinearSourceList &fit_sources,

        ///The list of sources dropped from the fit (need to be modified to
        ///set the amplitudes correctly for the best fit).
        LinearSourceList &dropped_sources,

        ///The PSF map to use in order to derive initial flux estimates.
        const PSF::PiecewiseBicubicMap &psf_map,

        ///The vertical grid boundaries defining the PSF cells.
        const std::vector<double> &x_grid,

        ///The horizontal grid boundaries defining the PSF cells.
        const std::vector<double> &y_grid,

        ///The sub-pixel sensitivity map.
        const Core::SubPixelMap &subpix_map,

        ///The maximum sum squared of amplitude changes allowed after a
        ///PSF/Amplitude fit iteration in order to consider the solution
        ///converged.
        double max_abs_amplitude_change,

        ///The maximum fraction the sum squared of the amplitude changes can
        ///be of the length of the vector of amplitudes after a PSF/Amplitude
        ///fit iteration in order to consider the solution converged.
        double max_rel_amplitude_change,

        ///Sources with reduced \f$\chi^2\f$ bigger than this are removed
        ///from the fit.
        double max_chi2,

        ///Pixels whose residuals after fitting is too large are not included
        ///in the fit. This parameter determines how large the residual
        ///cut-off is. See the help of --pix-rej command line option for
        ///details.
        double pixel_rejection,

        ///The minimum rate of convergence to require before giving up.
        double min_convergence_rate,

        ///The maximum number of iterations before giving up (use a negative
        ///value to disable this limit).
        int max_iterations,

        ///Smoothing penalty (see help for details).
        double smoothing_penalty,

        ///On output contains the best fit polynomial coefficients. Resized
        ///as necessary.
        Eigen::VectorXd &best_fit_poly_coef,

        ///The gain of the input image.
        double gain
    );

    ///Performs the combined amplitude fitting for a group of sources.
    void LIB_LOCAL fit_group(
        ///The group to fit.
        OverlapGroup<LinearSource, PSF::PiecewiseBicubic> &group,

        ///The sub-pixel sensitivity map.
        const Core::SubPixelMap &subpix_map,

        ///The best fit PSF map.
        const PSF::PiecewiseBicubicMap &psf_map
    );

    ///\brief Fits for the amplitudes of sources initially dropped from PSF
    ///fitting.
    ///
    ///Fills in the flux, its error estimate and \f$\chi^2\f$.
    void LIB_LOCAL fit_dropped_sources(
        ///The list of dropped sources.
        LinearSourceList &dropped_sources,

        ///the best fit polynomial coefficients produced by
        ///fit_piecewise_bicubic_psf()
        const Eigen::VectorXd &best_fit_poly_coef
    );

    ///\brief Generates a file containing all the information necessary in
    ///order to reproduce the best fit PSF for any image position.
    void LIB_LOCAL output_best_fit_psf(
        ///The best fit polynomial coefficients giving the expansion of all
        ///PSF parameters as a function of image position.
        const Eigen::VectorXd &best_fit_poly_coef,

        ///The vertical grid boundaries defining the PSF cells.
        const std::vector<double> &x_grid,

        ///The horizontal grid boundaries defining the PSF cells.
        const std::vector<double> &y_grid,

        ///The name to use for the file. Overwritten if already exists.
        const std::string &fname
    );

} //End FitPSF namespace.

#endif
