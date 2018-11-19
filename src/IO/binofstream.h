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
    private:
        ///The underlying stream to the open file.
        std::ofstream     m_outstream;

    public:
        ///\brief Make a binary output stream that will be bound to a file
        ///later using the open method inherited from std::ofstream
        binofstream() {}

        ///Make a stream that does binary output to the given file.
        binofstream(
            ///The name of the file to attach the stream to.
            const std::string& fname 
        )
        {
            m_outstream.open(fname.c_str());
        }

        ///Converable to the underlying output file stream.
        operator std::ofstream& ()
        {
            return m_outstream;
        }

        ///Converable to the underlying output file stream.
        operator const std::ofstream& () {
            return m_outstream;
        }

        ///Close the stream.
        void close()
        {
            m_outstream.close();
        }

    }; //End binofstream class.

} //End IO namespace.

#endif
