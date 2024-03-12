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

} //End Core namespace.
