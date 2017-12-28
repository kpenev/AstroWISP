#ifndef __SDK_PSF_FITTING_H
#define __SDK_PSF_FITTING_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SDKSourceBase.h"
#include "Common.h"
#include "Source.h"
#include "../PSF/EllipticalGaussian.h"
#include "../Background/Source.h"
#include "../Core/SourceID.h"
#include "../Core/Image.h"

namespace FitPSF {

    /**\brief A PSF fitting source class based on a LocalPolynomialPSF.
     *
     * \ingroup FitPSF
     */
    template<class SUBPIX_TYPE>
        class LIB_LOCAL SDKSource : public SDKSourceBase {
        private:
            ///\brief The integral of the normalized PSF over all pixels 
            ///belonging to the source
            std::list< std::valarray<double> > __psf_pixels;

            ///\brief An iterator over the PSF integrals and their 
            ///derivatives over pixels
            std::list< std::valarray<double> >::const_iterator __psf_int_iter;

            ///\brief Calculates the sub-pixel weighted integral of the given 
            ///PSF over a pixel.
            std::valarray<double> psf_integral(
                ///The left edge of the pixel relative to the source.
                double x,

                ///The lower edge of the pixel relative to the source.
                double y,

                ///The point spread function or pixel response function of
                ///the source (which is assumed depends in subpix_map).
                const PSF::EllipticalGaussian &psf,

                ///The sub-pixel sensitivity map. If the map resolution is
                ///0x0 the psf argument is assumed to be the PRF.
                const SUBPIX_TYPE &subpix_map,

                ///The relative precision to which the integra should be
                ///calculated.
                double precision
            );

        public:
            ///Identify the pixels which belong to a source
            template<class EIGEN_MATRIX>
                SDKSource(
                    ///The I.D. of the source for which to identify pixels.
                    const Core::SourceID &id,

                    ///The abscissa coordinate of the source center.
                    double x,

                    ///The oordinate coordinate of the source center.
                    double y,

                    ///The background under the pixel (in ADU)
                    const Background::Source &background,

                    ///The actual image we are deriving a PSF map for
                    const Core::Image<double> &observed_image,

                    ///See Source constructor.
                    const std::string &output_fname,

                    ///A two dimensional array which keeps track of what 
                    ///pixels of the input image are assigned to what source. 
                    ///On exit it is updated with the pixels belonging to the 
                    ///newly constructed source.
                    EIGEN_MATRIX &source_assignment,

                    ///The gain (electrons per ADU) in the observed image
                    double gain,

                    ///How much above the background a pixel needs to be in 
                    ///order to be allocated to this source (the alpha 
                    ///parameter in the description)
                    double alpha,

                    ///The id to assign to this source in the 
                    ///source_assignment array
                    int source_id,

                    ///If nonzero impose a circular aperture for the source
                    ///no larger than the given value (otherwise uses only
                    ///pixels inconsistent with the background at the
                    ///prescribed by alpha level). The size of the circular
                    ///aperture is the smallest size possible that
                    ///encapsulates all pixels that pass the alpha test.
                    double max_circular_aperture=0,

                    ///Should first order derivatives of the integral of the
                    ///PSF over a pixel be calculated.
                    bool calculate_first_deriv=false,

                    ///Should second order derivatives of the integral of the
                    ///PSF over a pixel be calculated.
                    bool calculate_second_deriv=false
                );

            using Source<PSF::EllipticalGaussian>::pixel_psf;

            ///\brief The integral of the normalized PSF over the current
            ///pixel and its derivatives.
            double pixel_psf(PSF::SDKDerivative deriv = PSF::NO_DERIV) const
            {return (*__psf_int_iter)[deriv];}

            ///\brief Recalculates the integrals of the PSF over the source
            ///pixels using the given PSF and sub-pixel map.
            void set_PSF(const PSF::EllipticalGaussian &psf,
                         const SUBPIX_TYPE &subpix_map);

            ///\brief Advance to the next pixel of the source returning true
            ///if it is not past the last pixel.
            virtual bool next_pixel()
            {
                return (Source<PSF::EllipticalGaussian>::next_pixel()
                        &&
                        ++__psf_int_iter != __psf_pixels.end());
            }

            ///Restarts the iteration over pixels from the beginning
            virtual void restart_pixel_iteration()
            {
                Source<PSF::EllipticalGaussian>::restart_pixel_iteration();
                __psf_int_iter = __psf_pixels.begin();
            }
        }; //End SDKSource class.

    typedef std::list< SDKSource<GSLSubPixType> >::iterator 
        GSLSourceIteratorType;

    ///\brief Calculates the sub-pixel weighted integral of the given PSF
    ///over the pixel with lower left corner at (x, y) relative to the source
    ///center. The presicion specified is relative.
    template<class SUBPIX_TYPE>
        std::valarray<double> SDKSource<SUBPIX_TYPE>::psf_integral(
            double x,
            double y,
            const PSF::EllipticalGaussian &psf,
            const SUBPIX_TYPE &subpix_map,
            double precision
        )
    {
        std::valarray<double> result(0.0, PSF::KK_DERIV + 1);
        if(subpix_map.x_resolution()==0) {
            assert(subpix_map.y_resolution()==0);
            psf.evaluate(x+0.5, y+0.5,
                         result,
                         calculate_first_deriv(),
                         calculate_second_deriv());
        } else {
            double dx=1.0/subpix_map.x_resolution(), 
                   dy=1.0/subpix_map.y_resolution(); 
            x+=0.5*dx; y+=0.5*dy;
            unsigned subpix_x_res=subpix_map.x_resolution();
            for(unsigned y_i=0; y_i<subpix_map.y_resolution(); y_i++)
                for(unsigned x_i=0; x_i<subpix_x_res; x_i++) {
                    psf.set_precision(precision, 0.0);
                    psf.integrate(
                            x+dx*x_i, y+dy*y_i,
                            dx, dy, 
                            0,
                            calculate_first_deriv(), calculate_second_deriv()
#ifdef DEBUG	
#ifdef SHOW_PSF_PIECES
                            , false, true
#endif
#endif
                    );
                    result[PSF::NO_DERIV] += (
                        psf.last_integrated(PSF::NO_DERIV)
                        *
                        subpix_map(x_i, y_i)
                    );
                    for(
                        int d = (calculate_first_deriv()
                                 ? PSF::S_DERIV
                                 : PSF::SS_DERIV);
                        d <= (calculate_second_deriv()
                              ? PSF::KK_DERIV
                              : PSF::K_DERIV);
                        ++d
                    ) {
                        result[d] += (
                            psf.last_integrated(
                                static_cast<PSF::SDKDerivative>(d)
                            )
                            *
                            subpix_map(x_i, y_i)
                        );
#ifdef DEBUG
                        assert(!std::isnan(result[d]));
#endif
                    }
                }
        }
        return result;
    }

    ///Identify the pixels which belong to the source centered at (x,y).
    template<class SUBPIX_TYPE>
        template<class EIGEN_MATRIX>
        SDKSource<SUBPIX_TYPE>::SDKSource(
            const Core::SourceID &id, double x, double y,

            ///The background under the pixel (in ADU)
            const Background::Source &background,

            ///The actual image we are deriving a PSF map for
            const Core::Image<double> &observed_image,

            ///See Source constructor.
            const std::string &output_fname,

            ///A two dimensional array which keeps track of what pixels of 
            ///the input image are assigned to what source. On exit it is 
            ///updated with the pixels belonging to the newly constructed 
            ///source.
            EIGEN_MATRIX &source_assignment,

            ///The gain (electrons per ADU) in the observed image
            double gain,

            ///How much above the background a pixel needs to be in order to
            ///be allocated to this source (the alpha parameter in the 
            ///description)
            double alpha,

            ///The id to assign to this source in the source_assignment array
            int source_id,

            ///If nonzero impose a circular aperture for the source no larger
            ///than the given value (otherwise uses only pixels inconsistent 
            ///with the background at the prescribed by alpha level). The 
            ///size of the circular aperture is the smallest size possible 
            ///that encapsulates all pixels that pass the alpha test.
            double max_circular_aperture,

            ///Should first order derivatives of the integral of the PSF
            ///over a pixel be calculated.
            bool calculate_first_deriv,

            ///Should second order derivatives of the integral of the PSF
            ///over a pixel be calculated.
            bool calculate_second_deriv
        ) :
            SDKSourceBase(id,
                          x,
                          y,
                          background,
                          observed_image,
                          source_assignment,
                          gain,
                          alpha,
                          source_id,
                          output_fname,
                          calculate_first_deriv,
                          calculate_second_deriv,
                          max_circular_aperture)
    {}

    ///\brief Recalculates the integrals of the PSF over the source pixels 
    ///using the given PSF and sub-pixel map.
    template<class SUBPIX_TYPE>
        void SDKSource<SUBPIX_TYPE>::set_PSF(
            const PSF::EllipticalGaussian &psf,
            const SUBPIX_TYPE &subpix_map
        )
        {
            Source<PSF::EllipticalGaussian>::restart_pixel_iteration();
            __psf_pixels.clear();
            do {
                __psf_pixels.push_back(pixel_psf(psf));
            } while(Source<PSF::EllipticalGaussian>::next_pixel());
            fit_flux();
        }

} //End FitPSF namespace.

#endif
