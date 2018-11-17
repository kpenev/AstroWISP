/**\file
 *
 * \brief Declares a class for working with single extension FITS images
 * with masks following fi conventions.
 *
 * \ingroup IO
 */

#ifndef __FITS_IMAGE_H
#define __FITS_IMAGE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Image.h"
#include "FitsHeader.h"
#include "fitsio.h"
#include "parse_hat_mask.h"
#include <typeinfo>

namespace IO {

    ///An interface for working with FITS images.
    template<typename DATA_TYPE>
        class LIB_PUBLIC FitsImage : public Core::Image<DATA_TYPE> {
        private:
            ///The value of each pixel. See values argument of Image::Image.
            DATA_TYPE *__pixel_values;

            ///Error estimateof each pixel. See errors argument of Image::Image.
            DATA_TYPE *__pixel_errors;

            ///The bad pixel mask. See mask argument of Image::Image.
            char *__mask;

            ///The FITS header
            FitsHeader __header;

            ///See filename().
            std::string __filename;

            ///Opens a fits file for reading.
            fitsfile *open_fits(const std::string &filename);

            ///Closes as previously opened fits file.
            void close_fits(
                fitsfile *fptr,
                const std::string &filename = "unknown fits file"
            );

            ///\brief Based on the fi function for parsing mask strings from
            ///the header.
            void parse_mask_string(
                ///The mask string to parse.
                const std::string &mask_string,

                ///The horizontal resolution of the image.
                long x_resolution,

                ///The vertical resolution of the image.
                long y_resolution
            ) {
                if(__mask) delete[] __mask;
                __mask = new char[x_resolution * y_resolution];
                memset(__mask, 0, x_resolution * y_resolution);
                if(!mask_string.empty())
                    parse_hat_mask(mask_string.c_str(),
                                   x_resolution,
                                   y_resolution,
                                   __mask);
            }

            ///\brief Read pixel data from an HDU in a FITS file.
            void input_image(
                ///The open fits file from which to read.
                fitsfile *fptr,

                ///The number of pixels to read from the image.
                long num_pixels,

                ///The memory location where to place the data.
                DATA_TYPE *destination,

                ///The HDU to read from. If zero, the currently selected HDU in
                ///fptr is read.
                unsigned read_hdu = 0
            );

        public:
            ///\brief Reads the data in the given fits image(s).
            ///
            ///If both expect_errors is false and error_image == "", the
            ///errors are estimated assuming Poisson statistics of the pixel
            ///fluxes.
            FitsImage(
                ///The filename of the image containing pixel fluxes as
                ///primary HDU.
                const std::string &filename = "",


                ///The HDU number to containing the required image. If,
                ///zero, the first non-trivial image HDU is used.
                unsigned values_hdu = 0,

                ///Should rounding real valued images to integer be allowed.
                ///Rounding is only performed if DATA_TYPE is integer and the
                ///input image contains floating point vaules.
                bool allow_rounding = false,

                ///The HDU number of the extension containing error estimates
                ///for the pixel values. If zero, no error estimate is read-in.
                unsigned errors_hdu = 0
            ) :
                __pixel_values(NULL),
                __mask(NULL)
            {
                if(filename.size())
                    open(filename, values_hdu, allow_rounding, errors_hdu);
            }

            ///Copies orig to *this.
            FitsImage(const FitsImage &orig)
                :
                    Core::Image<DATA_TYPE>(orig),
                    __header(orig.__header),
                    __filename(orig.__filename)
            {}

            ///\brief Open the given fits image (see constructor
            ///documentation for details).
            virtual void open(
                 ///The name of the file to open. Only used in error messages.
                const std::string &fits_filename,

                ///The HDU number to containing the pixel values.
                unsigned values_hdu = 0,

                ///Should rounding real valued images to integer be allowed.
                ///Rounding is only performed if DATA_TYPE is integer and the
                ///input image contains floating point vaules.
                bool allow_rounding = false,

                ///The HDU number to containing the estimated pixel errors. If
                ///zero, no error image is read.
                unsigned errors_hdu = 0
            );

            ///Returns immutable reference to the image header.
            virtual const FitsHeader &header() const {return __header;}

            ///Returns mutable reference to the image header.
            virtual FitsHeader &header() {return __header;}

            ///Copies RHS to *this.
            FitsImage &operator=(const FitsImage &rhs)
            {
                operator=(rhs);
                __header = rhs.__header;
                __filename = rhs.__filename;
                return *this;
            }

            ///The name of the FITS file to which this image is attached.
            const std::string &filename() const {return __filename;}

            ///Clean up any allocated arrays.
            ~FitsImage()
            {
                if(__pixel_values) delete[] __pixel_values;
                if(__mask) delete[] __mask;
            }
        }; //End FitsImage class.

    template<class DATA_TYPE>
        fitsfile *FitsImage<DATA_TYPE>::open_fits(
            const std::string &filename
        )
        {
            fitsfile *fptr;
            int fits_status=0;
            fits_open_file(&fptr, filename.c_str(), READONLY, &fits_status);
            if(fits_status) {
                std::ostringstream msg;
                msg << "Failed to open " << filename;
                throw Error::Fits(msg.str().c_str());
            }
            return fptr;
        }

    template<class DATA_TYPE>
        void FitsImage<DATA_TYPE>::close_fits(fitsfile *fptr,
                                              const std::string &filename)
        {
            int fits_status=0;
            fits_close_file(fptr, &fits_status);
            if(fits_status) {
                std::ostringstream msg;
                msg << "Failed to close " << filename;
                throw Error::Fits(msg.str().c_str());
            }
        }

    template<class DATA_TYPE>
        void FitsImage<DATA_TYPE>::input_image(fitsfile *fptr,
                                               long num_pixels,
                                               DATA_TYPE *destination,
                                               unsigned read_hdu)
        {
            long start_pixel[2]={1, 1};
            int data_type, has_nan, fits_status = 0;
            void *undefined_value=NULL;
            DATA_TYPE nan;
            if(std::numeric_limits<DATA_TYPE>::has_quiet_NaN) {
                nan = std::numeric_limits<DATA_TYPE>::quiet_NaN();
                undefined_value = &nan;
            }
            if(typeid(DATA_TYPE) == typeid(unsigned char))
                data_type = TBYTE;
            if(typeid(DATA_TYPE) == typeid(char))
                data_type = TSBYTE;
            else if(typeid(DATA_TYPE) == typeid(short))
                data_type = TSHORT;
            else if(typeid(DATA_TYPE) == typeid(unsigned short))
                data_type = TUSHORT;
            else if(typeid(DATA_TYPE) == typeid(int))
                data_type = TINT;
            else if(typeid(DATA_TYPE) == typeid(unsigned int))
                data_type = TUINT;
            else if(typeid(DATA_TYPE) == typeid(long))
                data_type = TLONG;
            else if(typeid(DATA_TYPE) == typeid(unsigned long))
                data_type = TULONG;
            else if(typeid(DATA_TYPE) == typeid(float))
                data_type = TFLOAT;
            else if(typeid(DATA_TYPE) == typeid(double))
                data_type = TDOUBLE;
            else throw Error::Type(
                "Unsupported DATA_TYPE for FitsImageData."
            );

            if(read_hdu > 0)
                fits_movabs_hdu(fptr, read_hdu, NULL, &fits_status);
            if(fits_status) {
                std::ostringstream msg;
                msg << "Failed to move to HDU #" << read_hdu
                    << " in FitsImageData::input_image.";
                throw Error::Fits(msg.str());
            }
            fits_read_pix(fptr,
                          data_type,
                          start_pixel,
                          num_pixels,
                          undefined_value,
                          destination,
                          &has_nan,
                          &fits_status);

            std::ostringstream error_message;
            error_message
                << "Failed to read image data in FitsImageData::read HDU #"
                << read_hdu
                << "!"
                << std::endl;
            if(fits_status) throw Error::Fits(error_message.str());

        }

    template<class DATA_TYPE>
        void FitsImage<DATA_TYPE>::open(
            const std::string &fits_filename,
            unsigned values_hdu,
            bool allow_rounding,
            unsigned errors_hdu
        )
        {
            __filename = fits_filename;
            std::cerr << "Set image name to " << __filename << std::endl;

            fitsfile *fptr = open_fits(fits_filename);

            int dimensions, fits_status = 0, bitpix, hdu_type = IMAGE_HDU;
            long naxes[2];
            for(
                ;
                fits_status == 0;
                fits_movrel_hdu(fptr, 1, &hdu_type, &fits_status)
            ) {
                fits_get_img_param(fptr,
                                   2,
                                   &bitpix,
                                   &dimensions,
                                   naxes,
                                   &fits_status);
                if(dimensions != 0 && hdu_type == IMAGE_HDU) break;
            }

#ifndef NDEBUG
            std::cerr << "Found primary image with resolution "
                      << naxes[0] << "x" << naxes[1]
                      << std::endl;
#endif

            assert(fits_status == 0);

            __header.read(fptr);

            parse_mask_string(__header["MASKINFO"],
                              naxes[0],
                              naxes[1]);


            if(fits_status)
                throw Error::Fits("Failed to read image parameters in "
                                  "FitsImageData::read.");
            if(dimensions!=2)
                throw Error::Fits(
                    "Only 2D image are supported at this time."
                );

            if(
                std::numeric_limits<DATA_TYPE>::is_integer
                &&
                bitpix < 0
                &&
                !allow_rounding
            )
                throw Error::Fits("Non-rounding integer FitsImage attached "
                                  "to a real valued image.");

            if(__pixel_values) delete[] __pixel_values;
            __pixel_values = new DATA_TYPE[naxes[0] * naxes[1]];

            input_image(fptr, naxes[0] * naxes[1], __pixel_values, values_hdu);

            if(errors_hdu > 0) {
                __pixel_errors = new DATA_TYPE[naxes[0] * naxes[1]];
                input_image(fptr,
                            naxes[0] * naxes[1],
                            __pixel_errors,
                            errors_hdu);
            } else
                __pixel_errors = NULL;

            close_fits(fptr, fits_filename);
            Core::Image<DATA_TYPE>::wrap(__pixel_values,
                                         __mask,
                                         naxes[0],
                                         naxes[1],
                                         __pixel_errors);
        }

} //End IO namespace.
#endif
