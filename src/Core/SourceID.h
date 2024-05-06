#ifndef __SOURCE_ID_H
#define __SOURCE_ID_H
/**\file
 *
 * \brief The declaration of the SourceID class.
 */

#include "../Core/SharedLibraryExportMacros.h"
#include "Error.h"
#include <string>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

namespace Core {

    ///The two components of the HAT i.d. (primary filed and source number)
    class LIB_PUBLIC SourceID {
    private:
        ///The string representation of the ID.
        std::string __id;

        ///\brief Is this an ID of the form HAT-FFF-SSSSSSS, which gets
        ///parsed into a source and a field.
        bool __is_hatid;

        ///The field entry (defined only if this is a HAT-ID).
        unsigned __field,

                 ///The source entry (defined only if this is a HAT-ID).
                 __source;

        ///Returns HAT field and source index, if ID is in HAT format:
        ///HAT-FFF-SSSSSSS (where F is the field, and S is the source index).
        ///F and S > 1 for valid IDs.
        bool parse_hatid();
    public:

        ///Construct out of the given string.
        SourceID( const std::string &id, bool assume_non_hat=false) :
            __id( id ),
            __is_hatid(!assume_non_hat && parse_hatid())
        {}

        ///Create a source ID from the given field and source number.
        SourceID(
            ///The HAT field the source belongs to
            unsigned field,

            ///The source number within the HAT field.
            unsigned source
        ) :
            __field(field),
            __source(source)
        {
            std::ostringstream temp;
            temp.fill('0');
            temp << "HAT-"
                 << std::setw(3) << field
                 << std::setw(7) << source;
            __id=temp.str();
        }

        ///Create a non-existent source
        SourceID() {}

        ///Copy constructor.
        SourceID(const SourceID& id) :
            __id(id.__id),
            __is_hatid(id.__is_hatid),
            __field(id.__field),
            __source(id.__source)
        {}

        ///Convert to string.
        virtual operator const std::string& () const {return __id;}

        ///Copy RHS to this.
        virtual SourceID &operator=(const SourceID &rhs)
        {
            __is_hatid = rhs.__is_hatid;
            __id = rhs.__id;
            __field = rhs.__field;
            __source = rhs.__source;
            return *this;
        }

        ///Copy RHS to this.
        virtual SourceID &operator=(const std::string &rhs)
        {__id=rhs; __is_hatid=parse_hatid(); return *this;}

        ///The string representation of the ID.
        virtual const std::string& str() const {return __id;}

        ///Self explanatory.
        virtual bool operator!= (const SourceID& rhs) const {
            return __id != rhs.__id;
        }

        ///Self explanatory.
        virtual bool operator== (const SourceID& rhs) const {
            return __id == rhs.__id;
        }

        ///Is this a HAT-ID?
        virtual bool is_hatid() const {return __is_hatid;}

        ///\brief The field corresponding to the given source.
        ///
        ///An exception is thrown if not a HAT-ID.
        virtual unsigned field() const
        {
            if(!__is_hatid)
                throw Error::Type("Asking for field of a non-HAT ID!");
            return __field;
        }

        ///\brief The source corresponding to the given source.
        ///
        ///An exception is thrown if not a HAT-ID.
        virtual unsigned source() const
        {
            if(!__is_hatid)
                throw Error::Type("Asking for source of a non-HAT ID!");
            return __source;
        }

        LIB_PUBLIC_IMPL friend std::istream& operator>>( std::istream& is, SourceID& id );


        LIB_PUBLIC_IMPL friend std::ostream& operator<<(std::ostream& os,
                                        const SourceID& id);
    }; //End SourceID class.

    ///Output the ID to the given stream.
    LIB_PUBLIC std::ostream& operator<<(
        ///The stream to write to.
        std::ostream& os,

        ///The source ID to write.
        const SourceID& id
    );

    ///Read the ID from the given stream.
    LIB_PUBLIC std::istream& operator>>(
        ///The stream to read from.
        std::istream& is,

        ///The variable to set to the parsed ID.
        SourceID& id
    );

} //End Core namespace.

#endif
