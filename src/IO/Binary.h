/***************************************************************************/
/* binary-io.h                                                             */
/***************************************************************************/

#ifndef	__BINARY_IO_H_INCLUDED
#define	__BINARY_IO_H_INCLUDED 1

#include "../Core/SharedLibraryExportMacros.h"
#include <stdio.h>

#define INT_ID 0
#define DOUBLE_ID 1

namespace IO {

    ///\brief Packs the given integer values to an array of characters, 
    ///along with all the information required to read that array from a 
    ///file.
    LIB_PUBLIC int int_to_binary(int *values,
                                 size_t val_count,
                                 unsigned long int *num_bytes,
                                 char **packed);

    /* Packs (up to the given precision) the given double values to an array of
     * characters, along with all the information required to read that array
     * from a file. */
    LIB_PUBLIC int double_to_binary(double *values,
                                    size_t val_count,
                                    double precision,
                                    int has_nan,
                                    unsigned long int *num_bytes,
                                    char **packed);

    /*Adds the given array of values to the file. If outputting doubles an extra
     * argument must be provided giving the desired precision. */
    LIB_PUBLIC int bin_output(FILE *outfile,
                              void *values,
                              size_t val_count,
                              int type_id, 
                              ...);

    /* Reads an array of values from the file and sets type_id according to the 
     * type of data read. If indices is not null, it should contain a sorted list
     * of numbers of length val_count, which select which values to return. If
     * indices is NULL all values are read and val_count is set to their number.
     * */
    LIB_PUBLIC void *bin_input(FILE *infile,
                               int *type_id,
                               unsigned long *val_count,
                               const unsigned long *indices);

    /* Advances the file position of the given stream to immediately after the
     * next column. */
    LIB_LOCAL int bin_input_skip(FILE *infile);

} //End IO namespace.

#endif
