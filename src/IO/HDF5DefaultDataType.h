/**\file
 * 
 * \brief Defines default datatypes to use for HDF5 files.
 *
 * \ingroup IO
 */

#ifndef __DEFAULT_HDF5_DATA_TYPE
#define __DEFAULT_HDF5_DATA_TYPE

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/Error.h"
#include <H5Cpp.h>
#include <vector>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <cmath>

namespace IO {

    template<typename SCALAR>
        class LIB_PUBLIC H5DefaultDataType {
        public:
            static const H5::PredType &native;
            inline static const H5::PredType &file(
                const std::vector<SCALAR> &value
            );
        };

    ///Less than comparison working on the absolute values.
    template<typename T>
        LIB_LOCAL bool less_than_abs(const T& a, const T &b)
        {return std::abs(a) < std::abs(b);}

    template<>
        inline const H5::PredType &H5DefaultDataType<double>::file(
            const std::vector<double> &
        )
        {
            return H5::PredType::IEEE_F32LE;
        }

    template<typename SCALAR>
        const H5::PredType &H5DefaultDataType<SCALAR>::file(
            const std::vector<SCALAR> &value
        )
        {
            if(!std::numeric_limits<SCALAR>::is_integer)
                throw Error::Type("Asking for default HDF5 file type of an "
                                  "unrecognized SCALAR type!");
            SCALAR max_value = std::abs(
                *std::max_element(value.begin(),
                                  value.end(),
                                  less_than_abs<SCALAR>)
            );
            bool has_sign = std::numeric_limits<SCALAR>::is_signed;
            if(max_value < 1 << (has_sign ? 7 : 8))
                return (has_sign
                        ? H5::PredType::STD_I8LE
                        : H5::PredType::STD_U8LE);
            else if(max_value < 1 << (has_sign ? 15 : 16)) 
                return (has_sign
                        ? H5::PredType::STD_I16LE
                        : H5::PredType::STD_U16LE);
            else if(max_value < 1 << (has_sign ? 31 : 32))
                return (has_sign
                        ? H5::PredType::STD_I32LE
                        : H5::PredType::STD_U32LE);
            else return (has_sign
                         ? H5::PredType::STD_I64LE
                         : H5::PredType::STD_U64LE);

        }

} //End IO namespace.

#endif
