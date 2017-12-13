/**\file
 *
 * \brief Defines an output stream for columns of numeric values as binary.
 *
 * \ingroup IO
 */

#ifndef __BINOSTREAM_H
#define __BINOSTREAM_H

#include "Binary.h"
#include <iostream>
#include <valarray>
#include <limits>

namespace IO {
    ///Performs binary output of column of integer or double values
    class binostream : virtual public std::ostream {
    private:
        ///\brief The next set of values to pack to the output stream, 
        ///already converted to properly scaled and offset integers.
        std::valarray<unsigned long> __converted;

        ///The offset to apply to the converted values before converting them 
        ///back to whatever type they were.
        long __offset;

        ///The maximum value of the converted column
        unsigned long __max_value;

        ///Whether to allow for NaN values
        bool __has_nan;

        ///The discretization used for real value conversion.
        double __discretization;

        ///\brief Converts a real-valued column to unsigned long integers 
        ///preserving the values to the currently set precision
        template<class REAL_TYPE>
            void convert_real(const std::valarray<REAL_TYPE> &column);

        ///\brief Converts an integer-valued column to unsigned long integers 
        ///losslessly
        template<class INT_TYPE>
            void convert_int(const std::valarray<INT_TYPE> &column);
    public:
        ///\brief Pack in binary form the given column, if VAL_TYPE is 
        ///integer the  packing is lossless, if it is real valued,
        ///the precision (number of digits after the decimal point to keep) 
        ///should be pre-set using the precision method inherited from 
        ///ostream.
        template<class VAL_TYPE>
            binostream &operator<<(const std::valarray<VAL_TYPE> &column);
    };

    template<class REAL_TYPE>
        void binostream::convert_real(const std::valarray<REAL_TYPE> &column)
        {
            double min_val = column.min(),
                   max_val = column.max(); 
            __discretization = std::pow(10.0,
                                        -static_cast<int>(precision()));
            __offset = static_cast<long>(
                floor(0.5 +  min_val / __discretization)
            );
            __max_value = (
                static_cast<unsigned long>(
                    floor(0.5 + max_val / __discretization)
                )
                -
                __offset
            );
            __has_nan = false;
            __converted.resize(column.size());
            for (size_t i = 0; i < column.size(); ++i)
                if (std::isnan(column[i]) || std::isinf(column[i])) {
                    __has_nan = 1;
                    __converted[i] = __max_value + 2;
                } else 
                    __converted[i] = (
                        static_cast<unsigned long>(
                            floor(0.5 + column[i] / __discretization)
                        )
                        -
                        __offset
                    );
        }

    template<class INT_TYPE>
        void binostream::convert_int(const std::valarray<INT_TYPE> &column)
        {
            __has_nan = false;
            __offset = column.min();
            __max_value = column.max() - __offset;
            __converted.resize(column.size());
            for (size_t i = 0; i < column.size(); ++i) 
                __converted[i] = column[i] - __offset;
        }

    template<class VAL_TYPE>
        binostream &binostream::operator<<(
            const std::valarray<VAL_TYPE> &column
        )
        {
            //	char type_id;
            char *packed;
            unsigned long num_bytes,
                          val_count = column.size();
            if(std::numeric_limits<VAL_TYPE>::is_integer) {
                /*		type_id=INT_ID; 
                        convert_int(column);*/
                std::valarray<int> int_column(column.size());
                for(size_t i = 0; i < column.size(); ++i)
                    int_column[i] = column[i];
                if(
                    int_to_binary(&(int_column[0]),
                                  val_count,
                                  &num_bytes,
                                  &packed)
                )
                    throw Error::Runtime(
                        "Failed to pack integer column to binary."
                    );
            } else {
                //{type_id=DOUBLE_ID; convert_real(column);}
                assert(sizeof(VAL_TYPE) == sizeof(double));
                int has_nan = 0;
                std::valarray<double> dbl_column(column.size());
                __discretization = std::pow(10.0,
                                            -static_cast<int>(precision()));
                for(size_t i = 0; i < column.size(); ++i) {
                    if(std::isnan(column[i]) || std::isinf(column[i])) 
                        has_nan = 1;
                    dbl_column[i] = column[i];
                }
                if(
                    double_to_binary(&(dbl_column[0]),
                                     val_count,
                                     __discretization,
                                     has_nan,
                                     &num_bytes,
                                     &packed)
                )
                    throw Error::Runtime(
                        "Failed to pack real value column to binary."
                    );
            }

            /*	pack(&__converted[0], val_count, type_id, __offset, __discretization, 
                __max_value, __has_nan, &num_bytes, &packed);*/
            write(packed, num_bytes);
            free(packed);
            return *this;
        }

} //End IO namespace.

#endif
