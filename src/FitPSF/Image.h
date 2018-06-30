/**\file
 *
 * \brief Defines a class describing an image of Pixel pixels.
 *
 * \ingroup FitPSF
 */

#ifndef __PSF_FIT_PIXEL_IMAGE_H
#define __PSF_FIT_PIXEL_IMAGE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Pixel.h"
#include "Source.h"
#include "Image.h"
#include "../Background/Source.h"
#include "../IO/FitsImage.h"
#include "../Core/NaN.h"

namespace FitPSF {

    /**\brief A class for managing the selection of pixels for PSF/PRF 
     * fitting.
     *
     * \ingroup FitPSF
     */
    template <class SOURCE_TYPE>
        class LIB_PUBLIC Image : public IO::FitsImage<double> {
        private:
            typedef std::vector< Pixel<SOURCE_TYPE>* > PixelVector;

            ///\brief Image of pointers to Pixels for tracking all pixels
            ///created for PSF fitting.
            PixelVector __fit_pixels;

            ///The gain to assume for the observed image.
            double __gain;

            ///\brief The value and variance in electrons and electrons^2
            ///respectively of a pixel.
            std::pair<double, double> value_variance(
                ///The x-coordinate of the pixel in the image.
                unsigned x,

                ///The y-coordinate of the pixel in the image.
                unsigned y
            ) const;

            ///Delete any pixels allocated by assign_to_source()
            void delete_allocated_pixels();

        public:
            ///\brief Create a fit pixel manager for the given image.
            Image(
                ///The filename of the image to fit.
                const std::string &filename = "",

                ///The HDU number to containing the required image. If,
                ///zero, the first non-trivial image HDU is used.
                unsigned hdu_number = 0,

                ///The gain to assume for the image.
                double gain = 1.0
            ) :
                IO::FitsImage<double>(filename, hdu_number),
                __fit_pixels(x_resolution() * y_resolution(), NULL),
                __gain(gain)
            {}

            ///\brief Create a fit pixel manager for the given image.
            Image(
                ///The filename of the image to fit.
                const std::string &filename,

                ///The HDU number to containing the image of pixel values. If,
                ///zero, the first non-trivial image HDU is used.
                unsigned values_hdu,

                ///The HDU containing the estimated standard deviation of each
                ///pixel.
                unsigned errors_hdu
            ) :
                IO::FitsImage<double>(filename, values_hdu, false, errors_hdu),
                __fit_pixels(x_resolution() * y_resolution(), NULL),
                __gain(1.0)
            {}

            ///\brief Create a fit pixel manager for the given image.
            ///
            ///The newly constructed objects holds an alias of the image
            ///data, so the input image should not be destructed before this
            ///one.
            Image(
                ///The image whose pixels will participate in PSF fitting.
                Core::Image<double> &observed_image,

                ///The gain to assume for the image.
                double gain
            ) :
                __fit_pixels(observed_image.x_resolution()
                             *
                             observed_image.y_resolution(),
                             NULL),
                __gain(gain)
            {
                Core::Image<double>::wrap(observed_image);
            }

            ///\brief Add the pixel at the given coordinates to the given
            ///source and return a pointer to the pixel.
            ///
            ///This is the only method that creates Pixel instances.
            Pixel<SOURCE_TYPE> *assign_to_source(
                ///The x-coordinate of the pixel to assign to the source.
                unsigned x,

                ///The y-coordinate of the pixel to assign to the source.
                unsigned y,

                ///The source to assign this pixel to.
                SOURCE_TYPE *source
            );

            ///Return the pixel at the given location.
            const Pixel<SOURCE_TYPE> *operator()(unsigned x,
                                                 unsigned y) const
            {
                assert(x < x_resolution());
                assert(y < y_resolution());
                assert(__fit_pixels[index(x, y)] != NULL);
                return __fit_pixels[index(x, y)];
            }

            ///\brief Return the signal to noise ratio with which a pixel
            ///sticks above a background.
            double background_excess(
                ///The x-coordinate of the pixel to return the signal to
                ///noise of.
                unsigned x,

                ///The y-coordinate of the pixel to return the signal to
                ///noise of.
                unsigned y,

                ///The background to assume under the pixel (in ADU)
                const Background::Source &background
            ) const;

            ///\brief Return the signal to noise ratio with which a pixel
            ///sticks above a background.
            double background_excess(
                ///The x-coordinate of the pixel to return the signal to
                ///noise of.
                unsigned x,

                ///The y-coordinate of the pixel to return the signal to
                ///noise of.
                unsigned y,

                ///The background value to assume under the pixel in
                ///electrons.
                double background_electrons,

                ///The background variance to assume untder the source in
                ///electrons^2.
                double background_electrons_variance
            ) const;

        //Intentionally hide IO::FitsImage::open, but disable clang warning.
#ifdef TOOLCHAIN_CLANG
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif
            ///See IO::FitsImage::open().
            virtual void open(
                ///The gain to assume for the image.
                double gain,

                ///See IO::FitsImage::open().
                const std::string &filename,

                ///See IO::FitsImage::open().
                unsigned hdu_number = 1,

                ///See IO::FitsImage::open().
                bool allow_rounding = false
            );

            ///\brief Wrap the given data in an image.
            ///
            ///See Core::Image::wrap() for description of the arguments.
            virtual void wrap(double *values,
                              char *mask,
                              unsigned long x_resolution,
                              unsigned long y_resolution,
                              double *errors)
            {
                IO::FitsImage<double>::wrap(values,
                                            mask,
                                            x_resolution,
                                            y_resolution,
                                            errors);
                __fit_pixels.resize(x_resolution * y_resolution);
                __gain = 1.0;
            }

            virtual void wrap(Core::Image<double> &image)
            {
                IO::FitsImage<double>::wrap(image);
                __fit_pixels.resize(image.x_resolution() * image.y_resolution());
                __gain = 1.0;
            }

            ///The gain of the image set at construction.
            double gain() const {return __gain;}

            ///Delete any pixels allocated through assign_to_source.
            virtual ~Image() {delete_allocated_pixels();}
        }; //End Image class.

#ifdef TOOLCHAIN_CLANG
    #pragma clang diagnostic pop
#endif

    template<class SOURCE_TYPE>
        Pixel<SOURCE_TYPE> *Image<SOURCE_TYPE>::assign_to_source(
            unsigned x,
            unsigned y,
            SOURCE_TYPE *source
        )
        {
            assert(x < x_resolution());
            assert(y < y_resolution());
#ifdef VERBOSE_DEBUG
            std::cerr << "Assigning Pixel("
                      << x
                      << ", "
                      << y
                      << ") to Source("
                      << source->x()
                      << ", "
                      << source->y()
                      << "): "
                      << source
                      << std::endl;
#endif
            typename PixelVector::iterator destination = (
                __fit_pixels.begin()
                +
                index(x, y)
            );

#ifdef VERBOSE_DEBUG
            std::cerr << "Pixel ";
#endif

            if(*destination == NULL) {
#ifdef VERBOSE_DEBUG
                std::cerr << "does not exist. ";
#endif
                *destination = new Pixel<SOURCE_TYPE>(x,
                                                      y,
                                                      value_variance(x, y),
                                                      photometry_flag(x, y),
                                                      source);
            } else {
#ifdef VERBOSE_DEBUG
                std::cerr << "exists. ";
#endif
                (*destination)->add_to_source(source);
            }

#ifdef VERBOSE_DEBUG
            std::cerr << "Now with "
                << (*destination)->sources().size()
                << " sources."
                << std::endl;
#endif

            return *destination;
        }

    template<class SOURCE_TYPE>
        std::pair<double, double> Image<SOURCE_TYPE>::value_variance(
            unsigned x,
            unsigned y
        ) const
        {
            double value = Core::Image<double>::operator()(x, y) * __gain;
#ifdef VERBOSE_DEBUG
            if(has_errors())
                std::cerr << "Variance from separate image = "
                          << std::pow(error(x, y) * __gain, 2)
                          << std::endl;
            else
                std::cerr << "Variance from poission statistics = "
                          << std::abs(value)
                          << std::endl;
#endif
            return std::pair<double, double>(
                value,
                (
                    has_errors()
                    ? std::pow(error(x, y) * __gain, 2)
                    : std::abs(value)
                )
            );
        }

    template<class SOURCE_TYPE>
        void Image<SOURCE_TYPE>::delete_allocated_pixels()
        {
#ifdef VERBOSE_DEBUG
            std::cerr << "Deleting all allocated FitPSF pixels." << std::endl;
#endif
            for(
                typename PixelVector::iterator p = __fit_pixels.begin();
                p != __fit_pixels.end();
                ++p
            )
                if(*p != NULL) delete *p;
        }

    template<class SOURCE_TYPE>
        double Image<SOURCE_TYPE>::background_excess(
            unsigned                    x,
            unsigned                    y,
            const Background::Source    &background
        ) const
        {
            std::pair<double, double> pixel_val_var = value_variance(x, y);
            return FitPSF::background_excess(pixel_val_var.first,
                                             pixel_val_var.second,
                                             background,
                                             __gain);
        }

    template<class SOURCE_TYPE>
        double Image<SOURCE_TYPE>::background_excess(
            unsigned  x,
            unsigned  y,
            double    background_electrons,
            double    background_electrons_variance
        ) const
        {
            std::pair<double, double> pixel_val_var = value_variance(x, y);
            return FitPSF::background_excess(pixel_val_var.first,
                                             pixel_val_var.second,
                                             background_electrons,
                                             background_electrons_variance);
        }

    template<class SOURCE_TYPE>
        void Image<SOURCE_TYPE>::open(double               gain,
                                      const std::string    &filename,
                                      unsigned             hdu_number,
                                      bool                 allow_rounding)
        {
            IO::FitsImage<double>::open(filename,
                                        hdu_number,
                                        allow_rounding);
#ifdef VERBOSE_DEBUG
            std::cerr << "Fit image: "
                      << IO::FitsImage<double>::filename()
                      << std::endl;
#endif
            delete_allocated_pixels();
            __fit_pixels.assign(x_resolution() * y_resolution(),
                                NULL);
            __gain = gain;
        }

} //End FitPSF namespace.

#endif
