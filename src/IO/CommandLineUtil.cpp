/**\file
 * \brief Implementation of tilities for command line parsing common to all 
 * tools.
 *
 * \ingroup IO
 */

#include "CommandLineUtil.h"

using namespace IO;

namespace IO {

    void verify_input_columns(const std::list<Phot::Columns> &columns,
                              bool need_PSF)
    {
        bool found_x = false,
             found_y = false,
             found_S = false,
             found_D = false, 
             found_K = false,
             found_id = false;
        for(
            std::list<Phot::Columns>::const_iterator ci=columns.begin();
            ci != columns.end();
            ++ci
        )
            switch(*ci) {
                case Phot::id : found_id = true; break;
                case Phot::x : found_x = true; break;
                case Phot::y : found_y = true; break;
                case Phot::S : found_S = true; break;
                case Phot::D : found_D = true; break;
                case Phot::K : found_K = true; break;
                default :;
            }
        if(!found_id)
            throw Error::CommandLine("Missing 'id' in input columns.");
        if(!found_x)
            throw Error::CommandLine("Missing 'x' in input columns.");
        if(!found_y)
            throw Error::CommandLine("Missing 'y' in input columns.");
        if(need_PSF) {
            if(!found_S)
                throw Error::CommandLine("Missing 'S' in input columns.");
            if(!found_D)
                throw Error::CommandLine("Missing 'D' in input columns.");
            if(!found_K)
                throw Error::CommandLine("Missing 'K' in input columns.");
        }
    }

} //End IO namespace.
