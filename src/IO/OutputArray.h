/**\file
 *
 * \brief Declares a class for writing array of values held in boost::any.
 *
 * \ingroup IO
 */

#ifndef __OUTPUT_ARRAY_H
#define __OUTPUT_ARRAY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "TranslateToAny.h"
#include "../Core/Typedefs.h"
#include "Eigen/Dense"
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <valarray>
#include <algorithm>

namespace IO {

    ///\brief Prepares to write an array of values from from boost::any.
    ///
    ///Each element of the array moust have UNIT_TYPE type, otherwise
    ///boost::bad_any_cast is thrown. The original data should have been
    ///either a std::vector<UNIT_TYPE> or an std::valarray<UNIT_TYPE>.
    template<typename UNIT_TYPE>
        class LIB_PUBLIC OutputArray {
        private:
            size_t __size;				///< The size of the array.
            const UNIT_TYPE *__data;	///< The first element.
            UNIT_TYPE *__allocated_data;///< Data allocated by this class.

            ///\brief Try reading the input data assuming it is in some type of
            ///container providing beging() and end() const_iterators
            template<class INPUT_ARRAY_TYPE>
                bool try_container_type(
                    ///The value to try parsing.
                    const boost::any &value
                );

            ///\brief Try reading the input data assuming it is in some type of
            ///container where entries are stored contigously in memory.
            template<class INPUT_ARRAY_TYPE>
                bool try_array_type(
                    ///The value to try parsing.
                    const boost::any &value
                );

        public:
            ///To be filled later using parse().
            OutputArray() : __allocated_data(NULL) {};

            ///Attempts all possible casts on the given value.
            OutputArray(
                ///The value(s) to initialize the array with.
                const boost::any &value
            ) : __allocated_data(NULL)
            {parse(value);}

            ///Parses the given value into this.
            void parse(
                ///The value to fill the array with.
                const boost::any &value
            );

            ///The number of elements in the array.
            const size_t &size() const {return __size;}

            ///\brief A pointer to the first element in the array, the rest
            ///are contiguous.
            const UNIT_TYPE *data() const {return __data;}

            ///Constant reference to an array element.
            const UNIT_TYPE &operator[](
                ///The index within the array to return.
                size_t index
            ) const
            {assert(index < __size); return __data[index];}

            ///\brief Compares two arrays element by element (empty arrays
            ///compare equal).
            bool operator==(
                ///The original array to copy.
                const OutputArray<UNIT_TYPE> &rhs
            );

            ~OutputArray() {if(__allocated_data) delete[] __allocated_data;}
        };

    template<typename UNIT_TYPE>
        template<class INPUT_ARRAY_TYPE>
        bool OutputArray<UNIT_TYPE>::try_container_type(const boost::any &value) ///This has different linkage need to figure out what exactly
        {
            try {
                const INPUT_ARRAY_TYPE &
                    input_array = TranslateToAny<INPUT_ARRAY_TYPE>().get_value(
                        value
                    );
                __allocated_data = new UNIT_TYPE[input_array.size()];
                std::copy(input_array.begin(),
                          input_array.end(),
                          __allocated_data);
                __size = input_array.size();
                __data = __allocated_data;
                return true;
            } catch(const boost::bad_any_cast &) {
                return false;
            }
        }

    template<typename UNIT_TYPE>
        template<class INPUT_ARRAY_TYPE>
        bool OutputArray<UNIT_TYPE>::try_array_type(const boost::any &value)
        {
            try {
                const INPUT_ARRAY_TYPE &
                    input_array = TranslateToAny<INPUT_ARRAY_TYPE>().get_value(
                        value
                    );
                __allocated_data = new UNIT_TYPE[input_array.size()];
                const UNIT_TYPE *start = &(input_array[0]);
                __size = input_array.size();
                __data = __allocated_data;

                std::copy(start, start + __size, __allocated_data);
                return true;
            } catch(const boost::bad_any_cast &) {
                return false;
            }
        }

    template<typename UNIT_TYPE>
        LIB_PUBLIC void OutputArray<UNIT_TYPE>::parse(const boost::any &value)
        {

            typedef Eigen::Matrix<UNIT_TYPE, Eigen::Dynamic, 1> vector_eigen;

            if(
                !(
                    try_container_type< std::vector<UNIT_TYPE> >(value)
                    ||
                    try_array_type< std::valarray<UNIT_TYPE> >(value)
                    ||
                    try_array_type< vector_eigen >(value)
                )
            )
                throw boost::bad_any_cast();
        }

    template<>
        void OutputArray<double>::parse(const boost::any &value);

    template<typename UNIT_TYPE>
        bool LIB_PUBLIC OutputArray<UNIT_TYPE>::operator==(const OutputArray &rhs)
        {
            for(size_t i = 0; i < __size; ++i)
                if(__data[i] != rhs.__data[i]) return false;
            return true;
        }

} //End IO namespace.

#endif
