/**\file
 *
 * \brief Stream class for compressed binary output of sources to file.
 *
 * \ingroup IO
 */

#ifndef __BINARY_IO_H
#define __BINARY_IO_H

#include "binostream.h"
#include <fstream>

namespace IO {

    ///A binary output stream connected to a file.
    class LIB_PUBLIC binofstream : public binostream {
    public:
        ///\brief Make a binary output stream that will be bound to a file
        ///later using the open method inherited from std::ofstream
        binofstream() {}

        ///Make a stream that does binary output to the given file.
        binofstream( const std::string& fname ) {
            m_outstream.open( fname.c_str() );
        }

        operator std::ofstream& () {
            return m_outstream;
        }

        operator const std::ofstream& () {
            return m_outstream;
        }

        void close() {
            m_outstream.close();
        }

    private:

        std::ofstream     m_outstream;
    }; //End binofstream class.

} //End IO namespace.

#endif
