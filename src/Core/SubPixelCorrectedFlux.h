/**\file
 *
 * \brief Defile a class that corrects for non-uniform subpixel sensitivity
 * and partial pixels in the aperture.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SUB_PIXEL_CORRECTED_FLUX_H
#define __SUB_PIXEL_CORRECTED_FLUX_H

#include "../Core/SharedLibraryExportMacros.h"
#include "FluxPair.h"
#include "../Core/Image.h"
#include "../Core/SubPixelMap.h"
#include "../PSF/PSF.h"
#include "../Core/Flux.h"
#include "../Core/NaN.h"

#include "Eigen/Dense"

#include <list>
#include <cmath>
#include <iostream>

namespace Core {

    ///\brief A class to do fixed aperture photometry, corrected for variable
    ///sensitivity of different parts of a pixel.
    ///
    ///\ingroup SubPixPhot
    template<class SUBPIX_TYPE>
    class LIB_PUBLIC SubPixelCorrectedFlux {
    private:
        ///The apertures over which to sum up the flux for each source.
        std::list<double> __apertures;

        ///The subpixel sensitivity map, normalized to average to 1.0
        const SUBPIX_TYPE *__subpixel_map;

        ///The actual pixel values of the fits image.
        const Core::Image<double> *__image;

        double
            ///Constant variance to add to flux dependent variance.
            __var_offset,

            ///See gain argument to constructor.
            __gain;

        ///\brief A set of sub-pixel structure independent matrices (one
        ///for each aperture), which can be used to apply a different
        ///sub-pixel sensitivity map to a previously extracted source flux.
        std::valarray<Eigen::MatrixXd> __subpix_scaling;

        ///The index within __subpix_scaling which needs to be filled next.
        std::valarray<int> __scaling_ind;

        ///\brief Whether to assume the presence of Poisson noise in pixel
        ///fluxes (disabled if pixel error estimates are available).
        bool __photon_noise;

        ///\brief Return the sum of the integrals of the PSF over each
        ///subpixel weighted and not weighted by the subpixel sensitivity.
        ///
        ///If __subpixel_map is NULL, the weighted flux it returns is 0.
        FluxPair pixel_fluxes(
            ///The left edge of the pixel relative to the PSF center.
            double x,

            ///The bottom edge of the pixel relative to the PSF center.
            double y,

            ///The PSF to assume.
            const PSF::PSF &psf,

            ///The radius of  the aperture.
            double r = Core::NaN,

            ///If fill_scaling is non-negative and the entry in __scaling_ind
            ///corresponding to the aperture with the given index, the row in
            ///the corresponding __subpix_scaling matrix indexed by that
            ///entry is filled with the integral of the PSF over individual
            ///sub-pixels. The proper entries can then be calculated by
            ///dividing by (observed counts)*(raw_flux).
            int fill_scaling = -1
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            ,
            bool reset_piece_id = false
#endif
#endif
        );

        ///\brief Return <0 if the pixel is entirely contained within the
        ///given aperture, >0 if it is entirely outside of it and 0 if it is
        ///partially inside and partially outside.
        int classify_pixel(double x, double y, double aperture) const;

        ///Estimate the variance based on the fluxes and measured counts.
        double variance(const FluxPair &fluxes,
                        double measured_counts) const;

        ///\brief Return true iff the pixel with the given index is inside
        ///the image.
        inline bool in_image(long x_i, long y_i) const;

        ///\brief Make sure all members necessary to derive sub-pixel scaling
        ///matrices are properly initialized.
        void init_subpix_scaling();

    public:
        ///\brief Create a flux measuring object that will use the given
        ///apertures.
        ///
        ///The subpixel sensitivity map and the image must be specified
        ///separately before use.
        SubPixelCorrectedFlux(
            double error_offset,
            const std::list<double> &apertures = std::list<double>(),
            double gain = 1.5,
            bool assume_no_photon_noise = false
        ) :
            __apertures(apertures),
            __subpixel_map(NULL),
            __image(NULL),
            __var_offset(error_offset * error_offset),
            __gain(gain),
            __subpix_scaling(apertures.size()),
            __scaling_ind(apertures.size()),
            __photon_noise(!assume_no_photon_noise)
        {}

        ///\brief Create a flux measuring object that will use the given
        ///apertures and subpixel sensitivity map.
        ///
        ///The image must be specified separately
        ///before use.
        SubPixelCorrectedFlux(
            const SUBPIX_TYPE &subpixel_map,
            double error_offset,
            const std::list<double> &apertures = std::list<double>(),
            double gain = 1.5,
            bool assume_no_photon_noise = false
        ) :
            __apertures(apertures),
            __subpixel_map(&subpixel_map),
            __image(NULL),
            __var_offset(error_offset * error_offset),
            __gain(gain),
            __subpix_scaling(apertures.size()),
            __scaling_ind(apertures.size()),
            __photon_noise(!assume_no_photon_noise)
        {}

        ///\brief Create a flux measuring object that will use the given
        ///apertures, subpixel sensitivity map and image.
        SubPixelCorrectedFlux(
            const Core::Image<double> &image,
            const SUBPIX_TYPE &subpixel_map,
            double error_offset,
            const std::list<double> &apertures = std::list<double>(),
            double gain = 1.5,
            bool assume_no_photon_noise = false
        ) :
            __apertures(apertures),
            __subpixel_map(&subpixel_map),
            __image(&(image)),
            __var_offset(error_offset*error_offset),
            __gain(gain),
            __subpix_scaling(apertures.size()),
            __scaling_ind(apertures.size()),
            __photon_noise(
                !(
                    assume_no_photon_noise
                    ||
                    image.has_errors()
                )
            )
        {}

        ///\brief Adds another aperture to the current list of apertures used
        ///to measure flux over.
        void add_aperture(double aperture) {__apertures.push_back(aperture);}

        ///\brief The number of apertures for which flux measurements will be
        ///performed.
        size_t number_apertures() {return __apertures.size();}

        ///Sets the subpixel sensitivity map to use.
        void set_subpixel_map(const SUBPIX_TYPE &subpix_map)
        {__subpixel_map = &subpix_map;}

        ///Specifies the image over which photometry will be performed.
        void set_image(const Core::Image<double> &image)
        {
            __image = &(image);
            if(image.has_errors()) __photon_noise = false;
        }

        ///\brief Actually measure the flux in the current set of apertures
        ///centered on the given coordinates.
        std::valarray<Core::Flux> operator()(
            ///The horizontal offset of the aperture center relative to the
            ///left edge of the image.
            double x,

            ///The vertical offset of the aperture center relative to the
            ///bottom edge of the image.
            double y,

            ///The PSF to assume.
            const PSF::PSF &psf,

            ///If NaN, the PSF should include a background and the returned
            //flux is not background subtracted. Otherwise, the PSF should
            //not include a background.
            double background = Core::NaN,

            ///An estimate for the error in the background.
            double background_error = Core::NaN,

            ///If true, a set of sub-pixel structure independent matrices
            ///(one per aperture) is prepared, which can be used with a
            ///vector of sub-pixel sensitivities to derive the proper
            ///corrected flux for the same source on the same  image with a
            ///different subpixel structure. These matrices can be accessed
            ///through the get_subpix_scaling method.
            bool prepare_subpix_scaling = false
        );

        ///\brief Return the last calculated subpixel independent scaling for
        ///the aperture with the given index.
        const Eigen::MatrixXd &get_subpix_scaling(size_t aperture_ind)
        {return __subpix_scaling[aperture_ind];}

        ///Set the gain to assume when deriving the flux.
        void set_gain(double gain) {__gain=gain;}

        ///Return the image on which the object is operating.
        const Core::Image<double> &image() {return *__image;}
    }; //End SubPixelCorrectedFlux class.

    ///\brief Apply a previously derived sub-pixel independent scaling to
    ///calculate the source flux corrected with the given sub-pixel
    ///sensitivity map.
    ///
    ///\return The sub-pixel sensitivity corrected flux.
    template<typename EIGEN_TP1, typename EIGEN_TP2>
        double apply_subpix_scaling(
            ///A matrix giving the sub-pixel independent scaling.
            const Eigen::MatrixBase<EIGEN_TP1> &scaling,

            ///The sub-pixel sensitivity map. The x sub-pixel index changes
            ///faster going down the vector.
            const Eigen::MatrixBase<EIGEN_TP2> &subpix_map
        )
        {return (scaling*subpix_map).array().inverse().sum();}

    ///See the two template version of apply_subpix_scaling() for details.
    template<typename EIGEN_TP>
        double apply_subpix_scaling(
            const Eigen::MatrixBase<EIGEN_TP> &scaling,
            const Core::SubPixelMap &subpix_map
        )
        {
            unsigned num_subpix=subpix_map.x_resolution()*subpix_map.y_resolution();
            Eigen::VectorXd subpix_vector(num_subpix);
            unsigned vec_i=0;
            for(unsigned y=0; y<subpix_map.y_resolution(); y++)
                for(unsigned x=0; x<subpix_map.x_resolution(); x++)
                    subpix_vector[vec_i++]=subpix_map(x,y);
            return apply_subpix_scaling(scaling, subpix_vector);
        }

    template<class SUBPIX_TYPE>
        FluxPair SubPixelCorrectedFlux<SUBPIX_TYPE>::pixel_fluxes(
            double x,
            double y,
            const PSF::PSF &psf,
            double r,
            int fill_scaling
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            ,
            bool reset_piece_id
#endif
#endif
        )
        {
            double dx = 1.0 / __subpixel_map->x_resolution(),
            dy = 1.0/__subpixel_map->y_resolution(),
            raw = 0.0,
            weighted = 0.0;
            x += 0.5 * dx;
            y += 0.5 * dy;
            unsigned subpix_x_res = __subpixel_map->x_resolution();
            for(unsigned y_i = 0; y_i < __subpixel_map->y_resolution(); y_i++)
                for(unsigned x_i = 0; x_i < subpix_x_res; x_i++) {
                    double
                        full_pix_int = psf.integrate(x + dx * x_i,
                                                     y + dy * y_i,
                                                     dx,
                                                     dy
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                                     ,
                                                     0,
                                                     reset_piece_id
#endif
#endif
                        ),
                        partial_pix_int = (
                            std::isnan(r)
                            ? full_pix_int
                            : psf.integrate(
                                x + dx * x_i,
                                y + dy * y_i,
                                dx,
                                dy,
                                r
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                ,
                                reset_piece_id
#endif
#endif
                            )
                        );
                    raw += partial_pix_int;
                    if(__subpixel_map)
                        weighted += full_pix_int * (*__subpixel_map)(x_i, y_i);
                    if(fill_scaling >= 0) {
                        if(__scaling_ind[fill_scaling] >= 0)
                            __subpix_scaling[fill_scaling](
                                __scaling_ind[fill_scaling],
                                x_i+y_i*subpix_x_res
                            ) = full_pix_int;
                    }
                }
            return FluxPair(raw, weighted);
        }

    template<class SUBPIX_TYPE>
        int SubPixelCorrectedFlux<SUBPIX_TYPE>::classify_pixel(
            double x,
            double y,
            double aperture
        ) const
        {
            double xlarge,
            xsmall,
            ylarge,
            ysmall,
            ap2 = aperture * aperture;
            if(x > 0.0) {xlarge = x + 1; xsmall = -x;}
            else if(x > -1.0) {xsmall = 0; xlarge = std::max(x + 1, -x);}
            else {xsmall = x + 1; xlarge = -x;}
            if(y > 0.0) {ylarge = y + 1; ysmall = -y;}
            else if(y > -1.0) {ysmall = 0; ylarge = std::max(y + 1, -y);}
            else {ysmall = y + 1; ylarge = -y;}
            if(xlarge * xlarge + ylarge * ylarge <= ap2) return -1;
            else if(xsmall * xsmall + ysmall * ysmall < ap2) return 0;
            else return 1;
        }

    template<class SUBPIX_TYPE>
        double SubPixelCorrectedFlux<SUBPIX_TYPE>::variance(
            const FluxPair &fluxes,
            double measured_counts_variance
        ) const
        {
            return (std::pow(fluxes.raw_flux() / fluxes.weighted_flux(), 2)
                    *
                    measured_counts_variance
                    +
                    __var_offset);
        }

    template<class SUBPIX_TYPE>
        inline bool SubPixelCorrectedFlux<SUBPIX_TYPE>::in_image(
            long x_i,
            long y_i
        ) const
        {
            return (
                x_i >= 0
                &&
                y_i >= 0
                &&
                x_i < static_cast<long>(__image->x_resolution())
                &&
                y_i<static_cast<long>(__image->y_resolution())
            );
        }

    template<class SUBPIX_TYPE>
        void SubPixelCorrectedFlux<SUBPIX_TYPE>::init_subpix_scaling()
        {
            if(__subpix_scaling.size() != __apertures.size()) {
                __subpix_scaling.resize(__apertures.size());
                __scaling_ind.resize(__apertures.size(), -1);
            } else __scaling_ind = -1;
            unsigned num_subpix = (__subpixel_map->x_resolution()
                                   *
                                   __subpixel_map->y_resolution());
            size_t i = 0;
            for(
                std::list<double>::const_iterator ap_i = __apertures.begin();
                ap_i != __apertures.end();
                ++ap_i
            )
                __subpix_scaling[i++].resize(
                    std::pow(std::ceil(2.0 * (*ap_i) + 1), 2),
                    num_subpix
                );
        }

    template<class SUBPIX_TYPE>
        std::valarray<Core::Flux>
        SubPixelCorrectedFlux<SUBPIX_TYPE>::operator()(
            double x,
            double y,
            const PSF::PSF &psf,
            double background,
            double background_error,
            bool prepare_subpix_scaling)
        {
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
            bool reset_piece_id = true;
#endif
#endif
            std::valarray<double>
                flux_values(0.0, __apertures.size()),
                flux_variances(0.0, __apertures.size()),
                total_psf_integral(0.0, __apertures.size()),
                undefined_pixel_integral(0.0, __apertures.size());
            std::valarray<Core::PhotometryFlag> quality_flags(
                Core::GOOD,
                __apertures.size()
            );
            size_t num_ap = __apertures.size();
            double largest_aperture = __apertures.back(),
                   dist_to_edge = std::min(x, y);
            dist_to_edge = std::min(dist_to_edge,
                                    __image->x_resolution() - x);
            dist_to_edge = std::min(dist_to_edge,
                                    __image->y_resolution() - y);

            if(dist_to_edge < largest_aperture) {
                std::list<double>::const_iterator ap_i = __apertures.end();
                --ap_i;
                for(
                    int flag_i = num_ap - 1;
                    *ap_i > dist_to_edge && flag_i >= 0;
                    --flag_i, --ap_i
                )
                    quality_flags[flag_i] = Core::BAD;
            }
            double min_var = 0;
            long max_dist = static_cast<long>(std::ceil(largest_aperture)),
                 x_center_pix = static_cast<long>(std::floor(x)),
                 y_center_pix = static_cast<long>(std::floor(y));
            if(std::isnan(background) || std::isnan(background_error))
                max_dist = -1;
            if(prepare_subpix_scaling) init_subpix_scaling();
            for(long dist = 0; dist <= max_dist; ++dist) {
                for(
                    long y_i = y_center_pix - dist;
                    y_i <= y_center_pix + dist;
                    ++y_i
                )
                    for(
                        long x_i = x_center_pix - dist;
                        x_i <= x_center_pix + dist;
                        x_i += (std::abs(y_i-y_center_pix)==dist
                                ? 1
                                : 2 * dist)
                    ) {
                        bool pix_in_image = in_image(x_i, y_i);
                        double measured_counts,
                               measured_counts_variance;
                        if(pix_in_image) {
                            measured_counts = (*__image)(x_i, y_i);
                            measured_counts_variance =
                                std::pow(background_error, 2);
                            if(__photon_noise)
                                measured_counts_variance +=
                                    measured_counts / __gain;
                            else if(__image->has_errors())
                                measured_counts_variance += std::pow(
                                    __image->error(x_i, y_i),
                                    2
                                );
                        } else continue;

                        int pixel_status = 1;
                        double pix_min_x = static_cast<double>(x_i) - x,
                               pix_min_y = static_cast<double>(y_i) - y;
                        FluxPair fluxes;
                        unsigned res_i = 0;
                        for(
                            std::list<double>::const_iterator
                                ap_i = __apertures.begin();
                            ap_i != __apertures.end();
                            ++ap_i
                        ) {
                            psf.set_precision(
                                0.01 * std::sqrt(measured_counts_variance)
                                /
                                std::abs(measured_counts)
                                ,
                                0.01 * std::sqrt(flux_variances[res_i])
                            );

                            if(pixel_status >= 0) {
                                pixel_status = classify_pixel(pix_min_x,
                                                              pix_min_y,
                                                              *ap_i);
                                int fill_scaling = -1;
                                if(
                                    prepare_subpix_scaling
                                    &&
                                    pixel_status <= 0
                                    &&
                                    pix_in_image
                                ) {
                                    fill_scaling = res_i;
                                    ++__scaling_ind[res_i];
                                }
                                Core::PhotometryFlag quality =
                                    __image->photometry_flag(x_i, y_i);
                                if(
                                    pixel_status <= 0
                                    &&
                                    quality > quality_flags[res_i]
                                )
                                    for(unsigned i = res_i; i < num_ap; ++i)
                                        quality_flags[i] = quality;
                                if(pixel_status < 0) {
                                    fluxes = pixel_fluxes(pix_min_x,
                                                          pix_min_y,
                                                          psf,
                                                          Core::NaN,
                                                          fill_scaling
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                                          , reset_piece_id
#endif
#endif
                                    );
                                    if(fill_scaling >= 0)
                                        for(
                                            unsigned i = res_i + 1;
                                            i < num_ap;
                                            ++i
                                        ) {
                                            __subpix_scaling[i].row(
                                                ++__scaling_ind[i]
                                            ) = __subpix_scaling[res_i].row(
                                                __scaling_ind[res_i]
                                            );
                                        }
                                } else if(pixel_status == 0) {
                                    fluxes = pixel_fluxes(pix_min_x,
                                                          pix_min_y,
                                                          psf,
                                                          *ap_i,
                                                          fill_scaling
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                                          , reset_piece_id
#endif
#endif
                                    );
                                }
#ifdef DEBUG
#ifdef SHOW_PSF_PIECES
                                reset_piece_id = false;
#endif
#endif
                            }
                            if(pixel_status <= 0) {
                                total_psf_integral[res_i] +=
                                    fluxes.raw_flux();
                                if(pix_in_image) {
                                    flux_values[res_i] += (
                                        measured_counts
                                        *
                                        fluxes.raw_flux()
                                        /
                                        fluxes.weighted_flux()
                                    );
                                    flux_variances[res_i] += variance(
                                        fluxes,
                                        measured_counts_variance
                                    );
                                    if(prepare_subpix_scaling)
                                        __subpix_scaling[res_i].row(
                                            __scaling_ind[res_i]
                                        ) /= (measured_counts
                                              *
                                              fluxes.raw_flux());
                                    if(
                                        res_i == 0
                                        ||
                                        flux_variances[res_i] < min_var
                                    )
                                        min_var = flux_variances[res_i];
                                } else undefined_pixel_integral[res_i] +=
                                    fluxes.raw_flux();
                            }
                            res_i++;
                        }
                    }
            }
            std::valarray<double> scalings = (
                total_psf_integral
                /
                (total_psf_integral - undefined_pixel_integral)
            );
            std::valarray<Core::Flux> result(num_ap);
            size_t i = 0;
            for(
                std::list<double>::const_iterator ap_i = __apertures.begin();
                ap_i != __apertures.end();
                ++ap_i
            ) {
                result[i].value() = (flux_values[i] * scalings[i]
                                     -
                                     M_PI * std::pow(*ap_i, 2) * background);
                if(prepare_subpix_scaling) {
                    __subpix_scaling[i] /= scalings[i];
                    __subpix_scaling[i].conservativeResize(
                        __scaling_ind[i] + 1,
                        Eigen::NoChange
                    );
                }
                result[i].error() = std::sqrt(flux_variances[i]
                                              *
                                              std::pow(scalings[i], 2));
                result[i].flag() = quality_flags[i];
                ++i;
            }
            return result;
        } //End SubPixelCorrectedFlux::operator() definition.

} //End SubPixPhot namespace.

#endif
