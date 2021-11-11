#include "Binary.h"
#include <limits>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#define MIN_ULONG_BYTES 4
#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))
#define MIN_MAX(values, val_count, min, max, i) \
	for (i=0; isnan((values)[i]) && i<(val_count-1); i++) {}\
	(min)=(values)[i];\
	(max)=(values)[i];\
	for (; i<(val_count); i++) {\
		if (! isnan((values)[i])) {\
			if ((values)[i]<(min)) (min)=(values)[i];\
			if ((values)[i]>(max)) (max)=(values)[i];\
		}\
	}\

#define MIN_MAX_NONAN(values, val_count, min, max, i) \
	i=0;\
	(min)=(values)[i];\
	(max)=(values)[i];\
	for (; i<(val_count); i++) {\
        if ((values)[i]<(min)) (min)=(values)[i];\
        if ((values)[i]>(max)) (max)=(values)[i];\
	}\

namespace IO {

    /*Converts the given array of doubles to appropriately scaled integers which
     * preserve the values up to the desired precision, also returns the offset
     * which was applied to the double values before rescaling and truncating,
     * and the maximum integer value in the resulting array (the minimum is zero
     * by definition).*/
    int double_to_int(double *values, size_t val_count, double precision,
            unsigned long int *converted, long *offset,
            unsigned long int *max, int *has_nan)
    {
        double min_val, max_val;
        size_t i;
        MIN_MAX(values, val_count, min_val, max_val, i);
        *offset=(long)(floor(0.5+min_val/precision));
        *max=(unsigned long int)(floor(0.5+max_val/precision))-*offset;
        *has_nan=0;
        for (i=0; i<val_count; i++)
            if (isnan(values[i]) || isinf(values[i])) {
                *has_nan=1;
                converted[i]=*max+2;
            } else converted[i]=(unsigned long int)(
                    floor(0.5+values[i]/precision))-*offset;
        return 0;
    }

    /*Returns the values that should be used to denote nan.*/
    unsigned long int nan_value(int num_bits)
    {
        return pow(2, num_bits)-1;
    }

    /*Adds a single value to dest with a given bit offset.*/
    void push(unsigned long val, int , int bit_offset, char *dest)
    {
        val<<=bit_offset;
        *(unsigned long *)(dest)|=val;
    }

    /*Extracts the first value occupying num_bits with the given bit offset from
     * source.*/
    unsigned long pop(int num_bits, int offset, char *source)
    {
        unsigned long val=*(unsigned long *)(source);
        int to_zero=sizeof(unsigned long)*8-num_bits;
        val=(val<<(to_zero-offset))>>to_zero;
        return val;
    }

    /*Packs the values to the given file  as tightly as possible.*/
    void pack(unsigned long *values, FILE *outfile, unsigned long val_count,
            unsigned long max_value, int has_nan)
    {
        int bit_offset=0;
        size_t dest=0, i;
        char num_bits;
        unsigned long int num_bytes, nan_repr;
        char *packed;

        if (max_value==0) num_bits=(has_nan ? 1 : 0);
        else num_bits=(short int)(floor(log((double)(max_value))/log(2.0)))+1;
        num_bytes=(size_t)(ceil((double)(num_bits*val_count)/8.0));
        if (has_nan) {
            nan_repr=nan_value(num_bits);
            max_value++;
        } else nan_repr=0;

        packed=(char*)(malloc(num_bytes));
        for(i=0; i<num_bytes; i++) packed[i]=0;
        for (i=0; i<val_count; i++) {
            if (values[i]>max_value)
                push(nan_repr, num_bits, bit_offset, &packed[dest]);
            else push(values[i], num_bits, bit_offset, &packed[dest]);
            dest+=(num_bits+bit_offset)/8;
            bit_offset=(bit_offset+num_bits)%8;
        }
        fwrite((char*)(&val_count), MIN_ULONG_BYTES, 1, outfile);
        fwrite((char*)(&num_bits), 1, 1, outfile);
        fwrite((char*)(&num_bytes), MIN_ULONG_BYTES, 1, outfile);
        fwrite(packed, num_bytes, 1, outfile);
        free(packed);
    }

    /* Returns the array of characters that should be written to a file for the
     * given truncated data.*/
    int bin_chararray(int type_id, int has_nan, long offset, double precision,
            unsigned long *truncated, unsigned long val_count,
            unsigned long max_value, char **packed, unsigned long *num_bytes)
    {
        int bit_offset=0;
        char num_bits;
        unsigned long int nan_repr, coded_offset;
        /* type id, precision, coded offset, number values, bits/value,
         * number bytes */
        size_t info_size=2 + (type_id==DOUBLE_ID ? sizeof(double) : 0) +
            3*MIN_ULONG_BYTES, val_dest=info_size, i;
        char *hdr_dest;


        if (max_value==0) num_bits=(has_nan ? 1 : 0);
        else num_bits=(short int)(floor(log((double)(max_value))/log(2.0)))+1;
        *num_bytes=(size_t)(ceil((double)(num_bits*val_count)/8.0));
        if (has_nan) {
            nan_repr=nan_value(num_bits);
            max_value++;
        } else nan_repr=0;

        *packed=(char*)(calloc(info_size+*num_bytes,1));
        if(*packed==NULL) return 1;

        /*add information section*/
        hdr_dest=*packed;
        (*packed)[0]=(char)(2*type_id+has_nan);
        hdr_dest++;
        if(type_id==DOUBLE_ID) {
            memcpy(hdr_dest, (void*)(&precision), sizeof(double));
            hdr_dest+=sizeof(double);
        }
        coded_offset=2*abs(offset)+(offset<0 ? 1 : 0);
        memcpy(hdr_dest, (void*)(&coded_offset), MIN_ULONG_BYTES);
        hdr_dest+=MIN_ULONG_BYTES;
        memcpy(hdr_dest, (void*)(&val_count), MIN_ULONG_BYTES);
        hdr_dest+=MIN_ULONG_BYTES;
        memcpy(hdr_dest, (void*)(&num_bits), 1);
        hdr_dest++;
        memcpy(hdr_dest, (void*)(num_bytes), MIN_ULONG_BYTES);
        hdr_dest+=MIN_ULONG_BYTES;

        /*add packed values*/
        for (i=0; i<val_count; i++) {
            if (truncated[i]>max_value)
                push(nan_repr, num_bits, bit_offset, &(*packed)[val_dest]);
            else push(truncated[i], num_bits, bit_offset, &(*packed)[val_dest]);
            val_dest+=(num_bits+bit_offset)/8;
            bit_offset=(bit_offset+num_bits)%8;
        }
        *num_bytes+=info_size;
        return 0;
    }

    /* Packs the given integer values to an array of characters, along with all
     * the information required to read that array from a file. */
    int int_to_binary(int *values, size_t val_count,
            unsigned long int *num_bytes, char **packed)
    {
        size_t i;
        int max_value, offset, result;
        unsigned long *truncated=
            (unsigned long *)(malloc(val_count*sizeof(unsigned long)));
        MIN_MAX_NONAN(values, val_count, offset, max_value, i);
        for (i=0; i<val_count; i++) truncated[i]=values[i]-offset;
        result=bin_chararray(INT_ID, 0, offset, 0, truncated, val_count,
                max_value-offset, packed, num_bytes);
        return result;
    }

    /* Packs (up to the given precision) the given double values to an array of
     * characters, along with all the information required to read that array
     * from a file. */
    int double_to_binary(double *values, size_t val_count, double precision,
            int has_nan, unsigned long int *num_bytes, char **packed)
    {
        unsigned long max_truncated;
        unsigned long *truncated=
            (unsigned long *)(malloc(val_count*sizeof(unsigned long)));
        long offset;
        if (double_to_int(values, val_count, precision, truncated, &offset,
                    &max_truncated, &has_nan)) return 1;
        if(bin_chararray(DOUBLE_ID, has_nan, offset, precision, truncated,
                    val_count, max_truncated, packed, num_bytes)) return 2;
        return 0;
    }

    /*Adds the given array of values to the file preserving values to
     * the given precision. The extra arguments are only used for floating point
     * values and they should be the desired precision and whether a special nan
     * value should be created.*/
    int bin_output(FILE *outfile, void *values, size_t val_count, int type_id,
            ...)
    {
        double precision;
        size_t i;
        long offset;
        unsigned long coded_offset;
        unsigned long int max_truncated;
        unsigned long int *truncated=(unsigned long *)(malloc(val_count*sizeof(unsigned long)));
        int has_nan, max_value;

        if (type_id==INT_ID) {
            has_nan=0;
            MIN_MAX_NONAN((int *)(values), val_count, offset, max_value, i);
            max_truncated=max_value-offset;
            for (i=0; i<val_count; i++)
                truncated[i]=((int *)(values))[i]-offset;
        } else {
            va_list prec_arg;
            va_start(prec_arg, type_id);
            precision=va_arg(prec_arg, double);
            has_nan=va_arg(prec_arg, int);
            va_end(prec_arg);
            if (double_to_int((double *)(values), val_count, precision,
                        truncated, &offset, &max_truncated, &has_nan)) return 1;
        }
        fputc((char)(2*type_id+has_nan), outfile);
        if (type_id==DOUBLE_ID)
            fwrite((char*)(&precision), sizeof(double), 1, outfile);
        coded_offset=2*abs(offset)+(offset<0 ? 1 : 0);
        if (!fwrite((char*)(&coded_offset), MIN_ULONG_BYTES, 1, outfile))
            return 1;
        pack(truncated, outfile, val_count, max_truncated, has_nan);
        free(truncated);
        return 0;
    }

    /* Reads an array of values from the file and sets type_id according to the
     * type of data read. If indices is not null, it should contain a sorted list
     * of numbers of length val_count, which select which values to return. If
     * indices is NULL all values are read and val_count is set to their number.
     * */
    void *bin_input(FILE *infile, int *type_id, unsigned long *val_count, const unsigned long *indices)
    {
        int has_nan;
        unsigned long num_bytes, source, i, result_ind;
        short int num_bits;
        double precision;
        long offset=0;
        unsigned long int nan_repr, int_val, array_size=0;
        int bit_offset=0;
        char *packed;
        void *result;

        *type_id=fgetc(infile);
        has_nan=*type_id%2;
        *type_id/=2;
        if (*type_id==DOUBLE_ID)
            if (!fread((char*)(&precision), sizeof(double), 1, infile))
                return NULL;
        int_val=0;
        if (!fread((char*)(&int_val), MIN_ULONG_BYTES, 1, infile)) return NULL;
        offset=(int_val%2 ? -1 : 1)*(long)(int_val)/2;
        if (!fread((char*)(&array_size), MIN_ULONG_BYTES, 1, infile))
            return NULL;
        if (indices==NULL) *val_count=array_size;
        num_bits=0;
        if (!fread((char*)(&num_bits), 1, 1, infile))
            return NULL;
        num_bytes=0;
        if (!fread((char*)(&num_bytes), MIN_ULONG_BYTES, 1, infile))
            return NULL;
        if (!num_bytes && num_bits) return NULL;
        nan_repr=nan_value(num_bits);
        packed=(char*)(malloc(num_bytes));
        if (num_bytes && !fread(packed, num_bytes, 1, infile)) return NULL;
        if (*type_id==INT_ID) result=malloc(*val_count*sizeof(int));
        else result=malloc(*val_count*sizeof(double));
        if (result==NULL) return NULL;
        source=0;
        result_ind=0;
        for (i=0; i<array_size; i++) {
            if (indices==NULL || i==indices[result_ind]) {
                if (!num_bits) int_val=0;
                else int_val=pop(num_bits, bit_offset, &packed[source]);
                if (int_val==nan_repr && has_nan)
                    ((double *)(result))[result_ind]=std::numeric_limits<double>::quiet_NaN();
                else if (*type_id==INT_ID)
                    ((int *)(result))[result_ind]=int_val+offset;
                else if (*type_id==DOUBLE_ID)
                    ((double *)(result))[result_ind]=precision*((long)int_val+
                                       offset);
                else {
                    free(packed);
                    free(result);
                    return NULL;
                }
                ++result_ind;
            }
            source+=(num_bits+bit_offset)/8;
            bit_offset=(bit_offset+num_bits)%8;
        }
        free(packed);
        return result;
    }

    int bin_input_skip(FILE *infile)
    {
        unsigned long num_bytes;
        char type_id;

        type_id=fgetc(infile);
        type_id/=2;
        if (type_id==DOUBLE_ID)
            if(fseek(infile, sizeof(double), SEEK_CUR)) return 1;
        if (fseek(infile, 2*MIN_ULONG_BYTES+1, SEEK_CUR)) return 2;
        num_bytes=0;
        if (!fread((char*)(&num_bytes), MIN_ULONG_BYTES, 1, infile)) return 3;
        if (fseek(infile, num_bytes, SEEK_CUR)) return 4;
        return 0;
    }

}//End IO namespace.
