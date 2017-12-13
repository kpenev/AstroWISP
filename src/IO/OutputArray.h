/**\file
 *
 * \brief Declares a class for writing array of values held in boost::any.
 *
 * \ingroup IO
 */

#ifndef __OUTPUT_ARRAY_H
#define __OUTPUT_ARRAY_H

#include "TranslateToAny.h"
#include "../Core/Typedefs.h"
#include "Eigen/Dense"
#include <H5Cpp.h>
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <valarray>

namespace IO {

    ///\brief Prepares to write an array of values from from boost::any.
    ///
    ///Each element of the array moust have UNIT_TYPE type, otherwise
    ///boost::bad_any_cast is thrown. The original data should have been 
    ///either a std::vector<UNIT_TYPE> or an std::valarray<UNIT_TYPE>.
    template<typename UNIT_TYPE>
        class OutputArray {
        private:
            hsize_t __size;				///< The size of the array.
            const UNIT_TYPE *__data;	///< The first element.
            UNIT_TYPE *__allocated_data;
        public:
            ///To be filled later using parse().
            OutputArray() : __allocated_data(NULL) {};

            ///Attempts all possible casts on the given value.
            OutputArray(const boost::any &value) : __allocated_data(NULL)
            {parse(value);}

            ///Parses the given value into this.
            void parse(const boost::any &value);

            ///The number of elements in the array.
            const hsize_t &size() const {return __size;}

            ///\brief A pointer to the first element in the array, the rest 
            ///are contiguous.
            const UNIT_TYPE *data() const {return __data;}

            ///Constant reference to an array element.
            const UNIT_TYPE &operator[](hsize_t index) const
            {assert(index < __size); return __data[index];}

            ///\brief Compares two arrays element by element (empty arrays 
            ///compare equal).
            bool operator==(const OutputArray<UNIT_TYPE> &rhs);

            ~OutputArray() {if(__allocated_data) delete[] __allocated_data;}
        };

    template<typename UNIT_TYPE>
        void OutputArray<UNIT_TYPE>::parse(const boost::any &value)
        {
            try {
                const std::vector<UNIT_TYPE>& vector =
                    TranslateToAny< std::vector<UNIT_TYPE> >().get_value(value);
                __size = vector.size();
                __data = &(vector[0]);
                return;
            } catch(boost::bad_any_cast) {}
            try {
                const std::valarray<UNIT_TYPE>& valarray =
                    TranslateToAny< std::valarray<UNIT_TYPE> >().get_value(value);
                __size = valarray.size();
                __data = &(valarray[0]);
                return;
            } catch(boost::bad_any_cast) {}
            typedef Eigen::Matrix<UNIT_TYPE, Eigen::Dynamic, 1> vector_eigen;
            const vector_eigen& 
                vector = TranslateToAny< vector_eigen >().get_value(value);
            __size = vector.size();
            __data = vector.data();
        }

    template<typename UNIT_TYPE>
        bool OutputArray<UNIT_TYPE>::operator==(const OutputArray &rhs)
        {
            for(hsize_t i = 0; i < __size; ++i)
                if(__data[i] != rhs.__data[i]) return false;
            return true;
        }

} //End IO namespace.

#endif
