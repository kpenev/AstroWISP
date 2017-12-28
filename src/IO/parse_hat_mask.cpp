/**\file
 *
 * \brief Define the parse_hat_mask function.
 *
 * \ingroup IO
 */

#include "parse_hat_mask.h"
#include <cstdio>

namespace IO {

    void parse_hat_mask(const char *mask_string,
                        long x_resolution,
                        long y_resolution,
                        char *mask)
    {
        const char FITS_MASK_MAX = 0x7F;
        const char FITS_MASK_DEF = 0x01;

        const char *next_term;
        static int xprev = 0,
                   yprev = 0,
                   use_diff = 0,
                   data = FITS_MASK_DEF;

        next_term = mask_string;
        while(*next_term) {
            int x, y, lx, ly;
            int num_scanned = std::sscanf(next_term,
                                          "%d,%d:%d,%d",
                                          &x, &y, &lx, &ly);
            if(num_scanned == 0)
                x = y = lx = ly = 0;
            else if(num_scanned == 1) {
                if(x > 0) use_diff = 1;
                else if(x<0) data=((-x) & FITS_MASK_MAX);
                else use_diff = 0;
                x = y = lx = ly = 0;
            } else if(num_scanned == 2)
                lx = ly = 1;
            else if(num_scanned == 3)
            {
                if(lx > 1)
                    ly = 1;
                else if(lx < -1) {
                    ly = -lx;
                    lx = 1;
                } else
                    lx = ly = 1;
            }

            if(lx > 0 && ly > 0)
            {
                if(use_diff) {
                    x += xprev;
                    y += yprev;
                }
                if(x < 0) {
                    lx += x;
                    x = 0;
                }
                if(y < 0) {
                    ly += y;
                    y = 0;
                }
                if( x + lx >= static_cast<int>(x_resolution) ) {
                    lx = x_resolution - x;
                }
                if( y + ly >= static_cast<int>(y_resolution) ) {
                    ly = y_resolution - y;
                }
                xprev = x;
                yprev = y;
                for( ; ly > 0 && lx > 0 ; y++, ly-- )
                {	
                    if( y < 0 || y >= static_cast<int>(y_resolution) )
                        continue;
                    char *to_update = mask + y * x_resolution + x;
                    for(int l = 0; l < lx; ++l) {
                        (*to_update) |= data;
                        ++to_update;
                    }
                }
            }
            while (*next_term && *next_term != ' ') ++next_term;
            while (*next_term == ' ') ++next_term;
        }
    }

}//End IO namespace.
