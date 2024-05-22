/**\file
 *
 * \brief Declare a base class for sources processed by SubPixTools.
 *
 * \ingroup Core
 */

#ifndef __SOURCE_H
#define __SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "SourceLocation.h"
#include "Flux.h"
#include "../Background/Source.h"
#include "../PSF/PSF.h"
#include <valarray>
#include <istream>
#include <limits>

namespace Core {

    ///A base class for all sources in AstroWISP.
    class LIB_PUBLIC Source : public SourceLocation {
    private:
        ///Is the source enabled (should it be processed?).
        bool __enabled;

        ///The flux estimate of the source for each aperture.
        std::valarray<Flux> __flux;

        ///An estimate of the backgound under the source.
        Background::Source __background;

        ///The point spread function applicable for this source.
        const PSF::PSF *__psf;

#ifdef DEBUG
        ///The time it took to process this source.
        double __processing_time;
#endif
    public:
        ///Create a source with the given id and position that will have its
        ///flux measured in the given number of apertures.
        Source(
            const SourceID& id,
            unsigned num_apertures = 1,
            double x = NaN,
            double y = NaN,
            const Background::Source &bg = Background::Source()
        ) :
            SourceLocation(id, x, y),
            __enabled(true),
            __flux(num_apertures),
            __background(bg),
            __psf(NULL) {}

        ///\brief Create a source with the given id and position that
        ///will have its flux measured in the given number of apertures, also
        ///specifying the point spread function that applies to this source.
        Source(
            unsigned num_apertures,
            const SourceID &id,
            double x,
            double y,
            const PSF::PSF &psf
        ) :
            SourceLocation(id, x, y),
            __flux(num_apertures),
            __psf(&psf)
        {}

        ///Copy orig to *this.
        Source(const Source &orig) :
            SourceLocation(orig),
            __flux(orig.__flux),
            __background(orig.__background),
            __psf(orig.__psf)
        {}

        ///\brief Sets the number of apertures in which this source's flux
        ///will be measured.
        void set_num_apertures(unsigned num_apertures)
        {__flux.resize(num_apertures);}

        ///\brief The estimate of the flux of this soruce in the aperture
        ///with the given index.
        const Flux &flux(unsigned aperture_index) const
        {return __flux[aperture_index];}

        ///\brief The estimate of the flux of this soruce in the aperture
        ///with the given index.
        Flux &flux(unsigned aperture_index) {return __flux[aperture_index];}

        ///\brief The estimate of the flux of this soruce in the aperture
        ///with the given index.
        const std::valarray<Flux> &flux() const {return __flux;}

        ///\brief The estimate of the flux of this soruce in the aperture
        ///with the given index.
        std::valarray<Flux> &flux() {return __flux;}

        ///The estimate of the background under this source.
        const Background::Source &background() const {return __background;}

        ///The estimate of the background under this source.
        Background::Source &background() {return __background;}

        ///The point spread function for this source
        const PSF::PSF &psf() const {return *__psf;}

        ///Change the point spread function for the source.
        void set_psf(const PSF::PSF &psf) {__psf=&psf;}

        ///Enables the source
        void enable() {__enabled=true;}

        ///Disables the source
        void disable() {__enabled=false;}

        ///Checks if the source is enabled
        bool is_enabled() const {return __enabled;}

#ifdef DEBUG
        ///The time it took to process this source.
        double processing_time() const {return __processing_time;}

        ///The time it took to process this source.
        double &processing_time() {return __processing_time;}
#endif

        ///\brief The S value for sources with elliptical Gaussian PSFs.
        ///
        ///An exception is thrown for non-elliptical Gaussian PSF sources.
        virtual double psf_s() const
        {throw Error::Type("Requesting the S value of a non-elliptical "
                           "Gaussian PSF source!");}

        ///\brief The D value for sources with elliptical Gaussian PSFs.
        ///
        ///An exception is thrown for non-elliptical Gaussian PSF sources.
        virtual double psf_d() const
        {throw Error::Type("Requesting the D value of a non-elliptical "
                           "Gaussian PSF source!");}

        ///\brief The K value for sources with elliptical Gaussian PSFs.
        ///
        ///An exception is thrown for non-elliptical Gaussian PSF sources.
        virtual double psf_k() const
        {throw Error::Type("Requesting the K value of a non-elliptical "
                           "Gaussian PSF source!");}

        ///\brief The best fit amplitude of the given source PSF in the
        ///image.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual double psf_amplitude() const
        {throw Error::Type("Requesting the PSF amplitude of a non-PSF "
                           "fitted source!");}

        ///\brief The best fit reduced \f$\chi^2\f$ of the given source in
        ///the image.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual double reduced_chi2() const
        {throw Error::Type("Requesting the reduced chi-squared of a non-PSF "
                           "fitted source!");}

        ///\brief Reference to the best fit reduced \f$\chi^2\f$ of the
        ///given source in the image.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual double &reduced_chi2()
        {throw Error::Type("Requesting a reference to the reduced "
                           "chi-squared of a non-PSF fitted source!");}

        ///\brief The signal to noise ratio of the given source in the
        ///image.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual double signal_to_noise() const
        {throw Error::Type("Requesting the signal to noise ratio of a "
                           "non-PSF fitted source!");}

        ///\brief Reference to the signal to noise ratio of the given source in
        ///the  image.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual double &signal_to_noise()
        {throw Error::Type("Requesting a reference to the signal to noise "
                           "ratio of a non-PSF fitted source!");}

        ///\brief The number of pixels assigned to the source.
        ///
        ///An exception is thrown for sources not generated by PSF fitting.
        virtual size_t flux_fit_pixel_count() const
        {throw Error::Type("Requesting the number of fitting pixels of a "
                           "non-PSF fitted source!");}

        ///Allow inheriting classes to define custom cleanup.
        virtual ~Source() {}

    }; //End Source class.

} //End Core namespace.

#endif
