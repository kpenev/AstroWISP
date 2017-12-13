/**\file
 *
 * \brief Declares a class for reading and working with FITS headers.
 *
 * \ingroup IO
 */

#ifndef __FITS_HEADER_H
#define __FITS_HEADER_H

#include "../Core/Error.h"
#include "fitsio.h"
#include <list>
#include <map>
#include <string>

namespace IO {

    ///\brief A structure representing the header of a fits file.
    ///
    ///\ingroup IO
    class FitsHeader {
    private:
        ///The keywords found in the header.
        std::list<std::string> __keywords;

        ///The value of the header keywords.
        std::map<std::string, std::string> __values;

        ///The comments of the header keywords
        std::map<std::string, std::string> __comments;
    public:
        ///Attach to the header in the given fits file
        FitsHeader(fitsfile *fptr=NULL)
        {if(fptr) read(fptr);}

        ///Copy orig to *this.
        FitsHeader(const FitsHeader &orig) :
            __keywords(orig.__keywords), __values(orig.__values) {}

        ///Attach to the header in the given fits file
        void read(fitsfile *fptr);

        ///Returns a list of the keywords in the header
        const std::list<std::string> &get_keywords() const
        {return __keywords;}

        ///Returns the value corresponding to the given keyword
        const std::string &get_value(const std::string &keyword)
        {return __values[keyword];}

        ///Returns the value corresponding to the given keyword
        const std::string &operator[](const std::string &keyword)
        {return __values[keyword];}

        ///Returns the comment to the given keyword.
        const std::string &get_comment(const std::string &keyword)
        {return __comments[keyword];}

        ///Copy rhs to *this.
        FitsHeader &operator=(const FitsHeader &rhs)
        {__keywords=rhs.__keywords; __values=rhs.__values; return *this;}
    };

} //End IO namespace.
#endif
