#include "Annulus.h"

namespace Background {

    Annulus::operator std::string() const
    {
        std::ostringstream result;
        result << __inner_radius << ":" << __width;
        return result.str();
    }

} //End Background namespace.
