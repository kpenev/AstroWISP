/**\file
 *
 * \brief Defines a class that allows placing inhomogeneous types data in
 * property trees.
 *
 * \ingroup IO
 */

#ifndef __TRANSLATE_TO_ANY_H
#define __TRANSLATE_TO_ANY_H

#include "../Core/SharedLibraryExportMacros.h"
#include <boost/any.hpp>
#include <string>

namespace IO {

    ///\brief Translator that works with boost::any for use with boost 
    ///property trees
    ///
    ///If a value to translate is passed by reference, it is copied (good for
    ///small values). If instead a pointer is passed, the pointer is stored, 
    ///but when the value is recovered the pointer is derereferenced. This 
    ///avoids copying the object, but means that the object must not be 
    ///destroyed (by e.g. going out of scope) before it is used.
    ///
    ///\ingroup IO
    template <typename DATA_TYPE>
        class LIB_LOCAL TranslateToAny
        {
        public:
            TranslateToAny() {}

            ///Convert the given value back to DATA_TYPE.
            const DATA_TYPE &get_value(
                ///The value to convert.
                const boost::any &value
            ) const
            {
                if(value.type()==typeid(DATA_TYPE))
                    return boost::any_cast< const DATA_TYPE& >(value);
                else return *boost::any_cast< const DATA_TYPE* >(value);
            }

            ///Construct a boost any-value containing a copy of the given value.
            boost::any put_value(
                ///The original value to copy.
                const DATA_TYPE& value
            ) const
            {return value;}

            ///\brief Construct a boost any-value containing a pointer to the 
            ///location pointed to by the input.
            boost::any put_value(
                ///The pointer to copy. 
                const DATA_TYPE* value
            ) const
            {return value;}
        };

    ///Translator for signed integers.
    const TranslateToAny<int> translate_int;

    ///Translator for unsigned integers.
    const TranslateToAny<unsigned> translate_unsigned;

    ///Translator for doubles.
    const TranslateToAny<double> translate_double;

    ///Translator for strings.
    const TranslateToAny<std::string> translate_string;

} //End IO namespace.

#endif
