#include "SourceID.h"

namespace Core {

    // returns true if a string contains digits only
    static bool is_digits( const char* str )
    {
        if ( !str || !*str ) {
            return false;
        }

        return str[ strspn( str, "0123456789" ) ] == '\0';
    }

    bool SourceID::parse_hatid()
    {
        // the HAT format must have a length of 15 characters
        if ( __id.length() != 15 ) {
            return false;
        }

        // buffers have the size of field (3) and source index length (7)
        // plus one, to accomodate the terminating byte
        char fieldbuf[4] = { '\0' };
        char sourcebuf[8] = { '\0' };
        if(
            std::sscanf( __id.c_str(), "HAT-%3c-%7c", fieldbuf, sourcebuf)
            !=
            2
        ) {
            return false;
        }

        // parse field ID
        if ( std::strlen( fieldbuf ) != 3 || !is_digits( fieldbuf ) ) {
            return false;
        }

        // parse source ID (sourcebuf has 7 chars length)
        if ( !is_digits( sourcebuf ) ) {
            return false;
        }

        // atoi cannot fail now
        __field = atoi( fieldbuf );
        __source = atoi( sourcebuf );
        return true;
    }

    std::istream& operator>>( std::istream& is, SourceID& id)
    {
        is >> id.__id;
        id.__is_hatid=id.parse_hatid();
        return is;
    }

    std::ostream& operator<<( std::ostream& os, const SourceID& id)
    {
        os << id.__id;
        return os;
    }

}
