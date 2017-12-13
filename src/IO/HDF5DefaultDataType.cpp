#include "HDF5DefaultDataType.h"

namespace IO {

    template<>
        const H5::PredType &H5DefaultDataType<char>::native =
        H5::PredType::NATIVE_CHAR;

    template<>
        const H5::PredType &H5DefaultDataType<unsigned char>::native =
        H5::PredType::NATIVE_UCHAR;

    template<>
        const H5::PredType &H5DefaultDataType<short>::native =
        H5::PredType::NATIVE_SHORT;

    template<>
        const H5::PredType &H5DefaultDataType<unsigned short>::native =
        H5::PredType::NATIVE_USHORT;

    template<>
        const H5::PredType &H5DefaultDataType<int>::native =
        H5::PredType::NATIVE_INT;

    template<>
        const H5::PredType &H5DefaultDataType<unsigned>::native =
        H5::PredType::NATIVE_UINT;

    template<>
        const H5::PredType &H5DefaultDataType<long>::native =
        H5::PredType::NATIVE_LONG;

    template<>
        const H5::PredType &H5DefaultDataType<unsigned long>::native =
        H5::PredType::NATIVE_ULONG;

    template<>
        const H5::PredType &H5DefaultDataType<double>::native =
        H5::PredType::NATIVE_DOUBLE;

} //End IO namespace.
