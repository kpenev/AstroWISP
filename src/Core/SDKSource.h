/**\file
 *
 * \brief Declare a class for sources with ellptical gaussian PSF.
 *
 * \ingroup Core
 */

#ifndef __SDK_SOURCE_H
#define __SDK_SOURCE_H

#include "Source.h"
#include "../PSF/EllipticalGaussian.h"

namespace Core {

    ///Represents a source with a PSF described by an elliptical gaussian.
    class LIB_PUBLIC SDKSource : public Source {
    private:
        double __s,///< The coefficient in front of (x^2+y^2) in the exponent
               __d,///< The coefficient in front of (x^2-y^2) in the exponent
               __k,///< The coefficient in front of (xy) in the exponent
               __amp,///< The amplitude of the PSF
               __bg;///< The background for the PSF only!

        ///The PSF of the source.
        PSF::EllipticalGaussian *__SDKpsf;
    public:
        ///\brief Create a source at the specified location with the specified
        ///PSF parameters, amplitude and background.
        SDKSource(
            const SourceID& id = SourceID(),
            unsigned num_apertures=1,
            double x=NaN,
            double y=NaN,
            double s=NaN,
            double d=NaN,
            double k=NaN,
            double amp=NaN,
            const Background::Source &bg=Background::Source(),
            double max_exp_coef=1
        ) :
            Source(id, num_apertures, x, y, bg), __SDKpsf(NULL)
        {set_psf(s, d, k, amp, bg.value(), max_exp_coef);}

        ///Copy a source.
        SDKSource(const SDKSource &orig) :
            Source(orig),
            __s(orig.__s), __d(orig.__d), __k(orig.__k), __amp(orig.__amp),
            __SDKpsf(new PSF::EllipticalGaussian(*orig.__SDKpsf))
        {Source::set_psf(*__SDKpsf);}

        ///The coefficient in front of (x^2+y^2) in the exponent
        double psf_s() const {return __s;}

        ///The coefficient in front of (x^2-y^2) in the exponent
        double psf_d() const {return __d;}

        ///The coefficient in front of (xy) in the exponent
        double psf_k() const {return __k;}

        ///The amplitude of the PSF
        double psf_amplitude() const {return __amp;}

        ///The backrgound under the PSF
        double psf_background() const {return __bg;}

        ///Changes the PSF of the source to have the given parameters.
        void set_psf(
            double s,
            double d,
            double k,
            double amp,
            double bg,
            double max_exp_coef = 1
        )
        {
            __s = s; __d = d; __k = k; __amp = amp, __bg = bg;
            if(__SDKpsf) {
                delete __SDKpsf;
            }
            __SDKpsf=new PSF::EllipticalGaussian(s,
                                               d,
                                               k,
                                               bg/amp,
                                               1e-5,
                                               0,
                                               max_exp_coef);
            Source::set_psf(*__SDKpsf);
        }

        ~SDKSource() {if(__SDKpsf) {delete __SDKpsf;}}
    }; //End SDKSource class.

} //End Core namespace.

///\brief Outputs all known information about an SDKSource object (i.e. the
///location and the PSF parameters).
std::ostream &operator<<(
    ///The stream to write to.
    std::ostream &os,

    ///The source to describe.
    const Core::SDKSource &src
);

#endif
