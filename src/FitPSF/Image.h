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
#include "../Core/Image.h"
#include "../Core/NaN.h"

namespace FitPSF {

    /**\brief A class for managing the selection of pixels for PSF/PRF
     * fitting.
     *
     * \ingroup FitPSF
     */
    template <class SOURCE_TYPE>
        class LIB_PUBLIC Image : public Core::Image<double> {
        private:
            typedef std::vector< Pixel<SOURCE_TYPE>* > PixelVector;

            ///\brief Image of pointers to Pixels for tracking all pixels
            ///created for PSF fitting.
            PixelVector __fit_pixels;

            ///The gain to assume for the observed image.
            double __gain;

            ///Unique identifier of this image.
            unsigned __id;

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
                ///Identifier to assign to this image.
                unsigned id = 0,

                ///The gain to assume for the image.
                double gain = 1.0
            ) :
                Core::Image<double>(),
                __gain(gain),
                __id(id)
            {
#ifdef VERBOSE_DEBUG
                std::cerr << "Created dummy FitPSF::Image instance at " << this
                          << std::endl;
#endif
            }

            unsigned id() const {return __id;}

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

#ifdef TOOLCHAIN_CLANG
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif TOOLCHAIN_GCC
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

            ///\brief Wrap the given data in an image.
            ///
            ///See Core::Image::wrap() for description of the arguments.
            virtual void wrap(
                ///The image values of the image to wrap.
                double *values,

                ///The pixel quality mask of the image to wrap.
                char *mask,

                ///The x resolution of the image.
                unsigned long x_resolution,

                ///The y resolution of the image.
                unsigned long y_resolution,

                ///Estimated errors of the image pixel values.
                double *errors,

                ///The ID to assign to the wrapped image.
                unsigned id
            )
            {
                Core::Image<double>::wrap(values,
                                            mask,
                                            x_resolution,
                                            y_resolution,
                                            errors);
                __fit_pixels.resize(x_resolution * y_resolution);
                __gain = 1.0;
                __id = id;
#ifdef VERBOSE_DEBUG
                std::cerr << "FitPSF::Image instance at " << this
                          << "with resolution "
                          << this->x_resolution() << "x" << this->y_resolution()
                          << " wrapped c-style data:"
                          << "\tvalues at: " << values
                          << "\terrors at: " << errors
                          << "\tmask at: " << (void*)mask
                          << std::endl;
#endif
            }

            ///Wrap the given image.
            virtual void wrap(Core::Image<double> &image)
            {
                Core::Image<double>::wrap(image);
                __fit_pixels.resize(image.x_resolution() * image.y_resolution());
                __gain = 1.0;
#ifdef VERBOSE_DEBUG
                std::cerr << "FitPSF::Image instance at " << this
                          << "with resolution "
                          << x_resolution() << "x" << y_resolution()
                          << " wrapped image at:" << &image
                          << std::endl;
#endif
            }

            ///The gain of the image set at construction.
            double gain() const {return __gain;}

            ///Delete any pixels allocated through assign_to_source.
            virtual ~Image() {
#ifdef VERBOSE_DEBUG
                std::cerr << "Destroying FitPSF::Image instance at " << this
                          << "with resolution "
                          << x_resolution() << "x" << y_resolution()
                          << std::endl;
#endif
                delete_allocated_pixels();
            }
        }; //End Image class.

#ifdef TOOLCHAIN_CLANG
    #pragma clang diagnostic pop
#elif TOOLCHAIN_GCC
    #pragma GCC diagnostic pop
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

} //End FitPSF namespace.

#endif
