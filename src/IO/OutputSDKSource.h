/**\file
 *
 * \brief The class of sources with Ellptical gaussian PSFs used for output.
 *
 * \ingroup IO
 */

#ifndef __OUTPUT_SDK_SOURCE_H
#define __OUTPUT_SDK_SOURCE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/SDKSource.h"

namespace IO {

    ///A class for SDK sources with all kinds of other attributes that can be
    ///set by the various fitting/photometering methods. For now additional
    ///information is:
    /// 1. Is the source enabled (for whatever)
    /// 2. The reduced chi2 assigned to this source by some fitting method
    /// 3. The number of pixels assigned to this source in an image
    class LIB_PUBLIC OutputSDKSource : public Core::SDKSource {
    private:
        ///The reduced \f$\chi^2\f$ for the source after PSF fitting.
        double __reduced_chi2,

               ///\brief The signal to noise in the pixels assigned to the 
               ///source (PSF fitting only).
               __signal_to_noise;

        ///The number of pixels assigned to the source (PSF fitting only).
        unsigned __num_pixels;

        std::string 
            ///See image_filename().
            __image_filename,

            ///See output_filename().
            __output_filename;

    public:
        OutputSDKSource(
                const std::string &image_fname,
                const std::string &output_fname,
                const Core::SourceID& id = Core::SourceID(),
                unsigned num_apertures = 1,
                double x = Core::NaN,
                double y = Core::NaN,
                double s = Core::NaN,
                double d = Core::NaN,
                double k = Core::NaN,
                double amp = Core::NaN,
                const Background::Source &bg = Background::Source(), 
                double reduced_chi2 = Core::NaN,
                unsigned num_pixels = 0,
                double signal_to_noise = Core::NaN
        ) :
            Core::SDKSource(id, num_apertures, x, y, s, d, k, amp, bg), 
            __reduced_chi2(reduced_chi2),
            __signal_to_noise(signal_to_noise), 
            __num_pixels(num_pixels),
            __image_filename(image_fname),
            __output_filename(output_fname)
        {}

        OutputSDKSource(const OutputSDKSource &orig) :
            Core::SDKSource(orig),
            __reduced_chi2(orig.__reduced_chi2),
            __signal_to_noise(orig.__signal_to_noise),
            __num_pixels(orig.__num_pixels),
            __image_filename(orig.__image_filename),
            __output_filename(orig.__output_filename)
        {}

        ///The reduced chi squared that is source had in the fit.
        double reduced_chi2() const {return __reduced_chi2;}

        ///The number of pixels assigned to this source during the fit
        unsigned pixel_count() const {return __num_pixels;}

        ///The total signal to noise of the source (all pixels combined).
        double signal_to_noise() const {return __signal_to_noise;}

        ///The filename of the image from which this source was extracted.
        const std::string &image_filename() const {return __image_filename;}

        ///The output filename where this source should be saved.
        const std::string &output_filename() const {return __output_filename;}

    }; //End OutputSDKSource class.

} //End IO namespace.

#endif
