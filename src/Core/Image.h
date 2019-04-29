/**\file
 *
 * \brief Declares a base class for images to be processed.
 *
 * \ingroup Core
 */

#ifndef __IMAGE_H
#define __IMAGE_H

#include "Typedefs.h"
#include "Error.h"
#include <cstddef>
#include <cstring>
#include <cassert>
#include <iostream>

namespace Core {

    const char MASK_OK              = 0x00; /* everything is ok... */
    const char MASK_CLEAR           = MASK_OK;  /* alias           */

    const char MASK_FAULT           = 0x01;/* faulty pixel         */
    const char MASK_HOT             = 0x02;/* hot (nonlinear) pixel*/
    const char MASK_COSMIC          = 0x04; /* hit by a cosmic ray */
    const char MASK_OUTER           = 0x08; /* outer pixel		  */
    const char MASK_OVERSATURATED   = 0x10; /* oversaturated       */
    const char MASK_LEAKED          = 0x20; /*leaked during readout*/
    const char MASK_SATURATED = (MASK_OVERSATURATED|MASK_LEAKED);
    const char MASK_INTERPOLATED    = 0x40; /*interpolated-not real*/

    const char MASK_BAD             = 0x7F; /* any error           */
    const char MASK_ALL             = MASK_BAD;

    const char MASK_NAN             = MASK_FAULT;

    ///Declare the minimum functionality expected from input images.
    template<typename DATA_TYPE>
        class LIB_LOCAL Image {
        private:
            ///\brief The pixel values in the image. See values argument of
            ///Image::Image()
            DATA_TYPE *__values;

            ///\brief Error estimates for __values. See errors argument of
            ///Image::Image()
            DATA_TYPE *__errors;

            ///\brief The quality flags of the image pixels. See mask
            ///argument of Image::Image()
            char *__mask;

            unsigned long
                ///The x resolutions of the image.
                __x_resolution,

                ///The y resolutions of the image.
                __y_resolution;

            ///\brief Is this image a wrap around another (and hence it should
            ///not free its memory).
            bool __wrapped;

        protected:
            ///\brief The index of the pixel at the given location within the
            ///1-D array.
            inline unsigned long index(unsigned long x,
                                       unsigned long y) const
            {
                assert(x < __x_resolution);
                assert(y < __y_resolution);
                return y * __x_resolution + x;
            }

            ///The pixel quality flag at the specified position.
            virtual char mask(unsigned long x, unsigned long y) const;

        public:
            ///See init_from_copy().
            Image(
                ///The pixel values in the image. The x coordinate changes
                ///faster than the y with successive entries. The image uses
                //a copy of this data.
                const DATA_TYPE *values = NULL,

                ///Quality flags for the image pixels in the same order as
                ///values. See comments for values argument.
                const char *mask = NULL,

                ///The number of pixels along the x direction in the image.
                unsigned long x_resolution = 0,

                ///The number of pixels along the y direction in the image.
                unsigned long y_resolution = 0,

                ///Error estimates for the image pixels.
                const DATA_TYPE *errors = NULL
            ) :
                __values(NULL),
                __errors(NULL),
                __mask(NULL),
                __x_resolution(0),
                __y_resolution(0),
                __wrapped(false)
            {
                if(values != NULL && x_resolution > 0 && y_resolution > 0)
                    init_from_copy(values,
                                   mask,
                                   x_resolution,
                                   y_resolution,
                                   errors);
            }

            ///\brief Create a new image containing a copy of the data of the
            ///given image.
            Image(
                ///The image to copy.
                const Image<DATA_TYPE> &image
            ) :
                __values(NULL),
                __errors(NULL),
                __mask(NULL),
                __x_resolution(0),
                __y_resolution(0),
                __wrapped(false)
            {
                if(image.__x_resolution > 0 && image.__y_resolution > 0) {
                    assert(image.__values);
                    init_from_copy(image.__values,
                                   image.__mask,
                                   image.__x_resolution,
                                   image.__y_resolution,
                                   image.__errors);
                }
            }

            ///Initialize this image with a copy of the given data.
            void init_from_copy(
                ///The pixel values in the image. The x coordinate changes
                ///faster than the y with successive entries. The image uses
                //a copy of this data.
                const DATA_TYPE *orig_values,

                ///Quality flags for the image pixels in the same order as
                ///values. See comments for values argument. If no bad pixel
                //mask should be used, set this to NULL.
                const char *orig_mask,

                ///The number of pixels along the x direction in the image.
                unsigned long orig_x_resolution,

                ///The number of pixels along the y direction in the image.
                unsigned long orig_y_resolution,

                ///The error estimates.
                const DATA_TYPE *orig_errors
            )
            {
                assert(__values == NULL);
                assert(__errors == NULL);
                assert(__mask == NULL);
                __values = new DATA_TYPE[orig_x_resolution
                                         *
                                         orig_y_resolution];
                __x_resolution = orig_x_resolution;
                __y_resolution = orig_y_resolution;
                std::memcpy(
                    __values,
                    orig_values,
                    orig_x_resolution * orig_y_resolution * sizeof(DATA_TYPE)
                );

                if(orig_mask) {
                    __mask = new char[orig_x_resolution * orig_y_resolution];
                    std::memcpy(__mask,
                                orig_mask,
                                orig_x_resolution * orig_y_resolution);
                } else
                    __mask = NULL;

                if(orig_errors) {
                    __errors = new DATA_TYPE[orig_x_resolution
                                             *
                                             orig_y_resolution];
                    std::memcpy(
                        __errors,
                        orig_errors,
                        (
                            orig_x_resolution
                            *
                            orig_y_resolution
                            *
                            sizeof(DATA_TYPE)
                        )
                    );
                } else
                    __errors = NULL;

#ifdef VERBOSE_DEBUG
                    std::cerr << "Copied to image at " << this
                              << " from" << std::endl
                              << "\tvalues = " << orig_values << std::endl
                              << "\terrors = " << orig_errors << std::endl
                              << "\tmask = " << (void*)orig_mask << std::endl
                              << "to" << std::endl
                              << "\t__values = " << __values << std::endl
                              << "\t__errors = " << __errors << std::endl
                              << "\t__mask = " << (void*)__mask << std::endl;
#endif

            }

            ///\brief Wrap the given data in an image.
            virtual void wrap(
                ///The pixel values in the image. The x coordinate changes
                ///faster than the y with successive entries. The data is
                ///used directly rather than copying, so it must not be
                ///destroyed while this object is still in use and my be
                ///modified through this object. To force copying the data,
                ///use the constructor with const arguments.
                DATA_TYPE *values,

                ///Quality flags for the image pixels in the same order as
                ///values. See comments for values argument.
                char *mask,

                ///The number of pixels along the x direction in the image.
                unsigned long x_resolution,

                ///The number of pixels along the y direction in the image.
                unsigned long y_resolution,

                ///Error estimates of the entries in values. If NULL, any query
                ///to the error results in an exception.
                DATA_TYPE *errors = NULL
            )
            {
                assert(__values == NULL);
                assert(__mask == NULL);
                assert(__errors == NULL);
                __values = values;
                __errors = errors;
                __mask = mask;
                __x_resolution = x_resolution;
                __y_resolution = y_resolution;
                __wrapped = true;
#ifdef VERBOSE_DEBUG
                std::cerr << "Wrapped image " << this
                          << " around:" << std::endl
                          << "\t__values = " << __values << std::endl
                          << "\t__errors = " << __errors << std::endl
                          << "\t__mask = " << (void*)__mask << std::endl;
#endif
            }

            ///Make this image an alias of the input image.
            virtual void wrap(
                ///The image to wrap
                Image<DATA_TYPE> &image
            )
            {
                wrap(image.__values,
                     image.__mask,
                     image.__x_resolution,
                     image.__y_resolution,
                     image.__errors);
            }

            ///\brief Was an error estimate provided for this image.
            bool has_errors() const {return __errors != NULL;}

            ///\brief An estimate of the error of the given pixel.
            DATA_TYPE &error(unsigned long x, unsigned long y)
            {
                assert(__errors != NULL);
                return __errors[index(x, y)];
            }

            ///\brief A place-holder for future development where each pixel
            ///may have an error estimate.
            DATA_TYPE error(
                ///The x coordinate of the pixel whose error to return.
                unsigned long x,

                ///The y coordinate of the pixel whose error to return.
                unsigned long y) const
            {
                assert(__errors != NULL);
                return __errors[index(x, y)];
            }

            ///The number of pixels in the x direction.
            unsigned long x_resolution() const {return __x_resolution;}

            ///The number of pixels in the y direction.
            unsigned long y_resolution() const {return __y_resolution;}

            ///Get the pixel value at the specified position.
            DATA_TYPE operator()(unsigned long x, unsigned long y) const
            {return __values[index(x, y)];}

            ///Reference to the pixel value at the specified position.
            DATA_TYPE &operator()(unsigned long x, unsigned long y)
            {
                return __values[index(x, y)];
            }

            ///True if the given pixel has no mask flags set.
            virtual bool good(unsigned long x, unsigned long y) const
            {return mask(x, y) == MASK_CLEAR;}

            ///True if the given pixel is flagged as faulty in the mask.
            virtual bool fault(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_FAULT;}

            ///True if the given pixel is flagged as hot in the mask.
            virtual bool hot(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_HOT;}

            ///True iff the given pixel is flagged as hit by a cosmic ray.
            virtual bool cosmic(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_COSMIC;}

            ///\brief True if the given pixel is flagged as coming from
            ///outside the image.
            ///TODO: this flag is never set
            virtual bool outside(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_OUTER;}

            ///True if the given pixel is flagged as oversaturated.
            virtual bool oversaturated(unsigned long x,
                                       unsigned long y) const
            {return mask(x, y) & MASK_OVERSATURATED;}

            ///True if the given pixel is flagged as being leaked to.
            virtual bool leaked(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_LEAKED;}

            ///True if the given pixel is oversaturated or leaked to.
            virtual bool saturated(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_SATURATED;}

            ///\brief True iff the given pixel is flagged having its value
            ///interpolated (rather than really read out).
            virtual bool interpolated(unsigned long x, unsigned long y) const
            {return mask(x,y) & MASK_INTERPOLATED;}

            ///True if any flag is set in the mask for the given pixel.
            virtual bool bad(unsigned long x, unsigned long y) const
            {return mask(x, y) & MASK_BAD;}

            ///\brief Return the photometry quality flag to assign to results
            ///involving the given pixel.
            virtual Core::PhotometryFlag photometry_flag(
                unsigned long x,
                unsigned long y
            ) const;

            virtual ~Image()
            {
                if(!__wrapped) {
#ifdef VERBOSE_DEBUG
                    std::cerr << "Deleting image at " << this
                              << " with" << std::endl
                              << "\t__values = " << __values << std::endl
                              << "\t__errors = " << __errors << std::endl
                              << "\t__mask = " << (void*)__mask << std::endl;
#endif
                    if(__values) delete[] __values;
                    if(__errors) delete[] __errors;
                    if(__mask) delete[] __mask;
                }
            }
        };

    template<typename DATA_TYPE>
        char Image<DATA_TYPE>::mask(unsigned long x, unsigned long y) const
        {
            if(__mask == NULL)
                return MASK_CLEAR;
            else if(
                x < 0 || y < 0 || x > x_resolution() || y > y_resolution()
            )
                return MASK_OUTER;
            else
                return __mask[index(x, y)];
        }

    template<typename DATA_TYPE>
        Core::PhotometryFlag Image<DATA_TYPE>::photometry_flag(
            unsigned long x,
            unsigned long y
        ) const
        {
            if(
                fault(x, y)
                ||
                hot(x, y)
                ||
                cosmic(x, y)
                ||
                interpolated(x, y)
                ||
                outside(x, y)
            )
                return Core::BAD;
            else if(saturated(x, y))
                return Core::SATURATED;
            else {
                assert(good(x, y));
                return Core::GOOD;
            }
        }
} //End Core namespace.
#endif
