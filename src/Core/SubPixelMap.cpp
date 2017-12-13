#include "SubPixelMap.h"

///Outputs the sensitivities of the subpixels as an array to the given stream.
std::ostream &operator<<(std::ostream &os,
                         const Core::SubPixelMap &subpix_map)
{
	std::ios_base::fmtflags orig_flags=os.flags();
	std::streamsize orig_precision=os.precision();
	os.precision(5);
	os.setf(std::ios::scientific);
	for(long y=subpix_map.y_resolution()-1; y>=0; y--) {
		for(unsigned long x=0; x<subpix_map.x_resolution()-1; x++)
			os << subpix_map(x,y) << ' ';
		os << subpix_map(subpix_map.x_resolution()-1,y) << std::endl;
	}
	os.flags(orig_flags);
	os.precision(orig_precision);
	return os;
}

namespace Core {

    ///Saves the sub-pixel map to a FITS image with the given name.
    void SubPixelMap::save_to_fits(const std::string &filename)
    {
        fitsfile *fptr;
        int fits_status = 0;
        fits_create_file(&fptr, filename.c_str(), &fits_status);
        if(fits_status)
            throw Error::Fits(("Failed to create " + filename).c_str());
        long naxes[] = {__x_res, __y_res};
        fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &fits_status);
        if(fits_status)
            throw Error::Fits(
                ("Failed to create image HDU in " + filename).c_str()
            );
        long first_pixel[] = {1, 1};
        fits_write_pix(fptr,
                       TDOUBLE,
                       first_pixel,
                       __x_res * __y_res,
                       &(__sensitivities[0]),
                       &fits_status);
        if(fits_status)
            throw Error::Fits(
                (
                    "Failed to write sub-pixel sensitivities to the first "
                    "HDU (Image) in " + filename
                ).c_str()
            );
        fits_close_file(fptr, &fits_status);
        if(fits_status)
            throw Error::Fits(("Failed to close " + filename).c_str());
    }

} //End Core namespace.
