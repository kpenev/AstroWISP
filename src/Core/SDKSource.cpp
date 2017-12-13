#include "SDKSource.h"

///Outputs all known information about an SDKSource object (i.e. the location
///and the PSF parameters).
std::ostream &operator<<(std::ostream &os, const Core::SDKSource &src)
{
	os << "Source at (" << src.x() << ", " << src.y() 
		<< ") with (S, D, K, Amp, Bg)=(" << src.psf_s() << ", " 
		<< src.psf_d() << ", " << src.psf_k() << ", " 
		<< src.psf_amplitude() << ", " << src.psf_background() << ")";
	return os;
}
