/**\file
 *
 * \brief Declare an input stream for sources class.
 *
 * \ingroup SubPixPhot
 */

#ifndef __SOURCE_INPUT_H
#define __SOURCE_INPUT_H

#include "../Core/PhotColumns.h"
#include "../Core/SDKSource.h"
#include "../IO/OutputSDKSource.h"
#include <iostream>
#include <list>

namespace SubPixPhot {

    ///A class that behaves as an input stream for sources
    class SourceInput {
    private:
        ///The stream from which input will be done
        std::istream *__instream;

        int __id_col,///< The column number of the source id
            __x_col, ///< The column number of the source x position
            __y_col, ///< The column number of the source y position
            __s_col, ///< The column number of the source S parameter
            __d_col, ///< The column number of the source D position
            __k_col, ///< The column number of the source K position
            __amp_col, ///< The column number of the source amplitude
            __bg_col, ///< The column number of the source background

            ///\brief The column number which specifies if the source is 
            ///enabled (for OutputSDKSource input only).
            __on_col, 

            __col_num;///< The number of columns in an input line

        ///\brief The minimum subdivision size to allow the read sources when
        ///they calculate intergrals of their PSFs
        double __max_exp_coef;

        bool __good, ///< Did the last operation complete successfully
             __eof; ///< Did we reach the end of input.

        ///\brief Set the values of the __*_col members according to the 
        ///contents of the argument.
        void set_columns(const std::list<Phot::Columns> &columns);

        ///\brief Read a source from the input stream to source and return 
        ///its enabled/disabled status if such a column was specified (the 
        ///result is undefined if it was not).
        template<class SOURCE_TYPE>
            bool read_source(SOURCE_TYPE &source,
                             double &s,
                             double &d,
                             double &k,
                             double &amp,
                             double &bg);
    public:
        ///\brief Construct an input stream that expects the input to contain 
        ///the given columns in the specified order.
        SourceInput(const std::list<Phot::Columns> &columns,
                    double max_exp_coef = 1) :
            __instream(NULL),
            __col_num(columns.size()),
            __max_exp_coef(max_exp_coef),
            __good(true),
            __eof(false) 
        {set_columns(columns);}

        ///Create a source input tied to the given stream.
        SourceInput(const std::list<Phot::Columns> &columns,
                    std::istream &instream,
                    double max_exp_coef = 1) :
            __instream(&instream),
            __col_num(columns.size()),
            __max_exp_coef(max_exp_coef),
            __good(instream),
            __eof(instream.eof())
        {set_columns(columns);}

        ///Read sources from the given stream
        void attach_to_stream(std::istream &instream)
        {__instream = &instream;}

        ///Read a source with an Elliptical Gaussian PSF.
        SourceInput &operator>>(Core::SDKSource &source);

        ///Read a source with an unspecified PSF.
        SourceInput &operator>>(Core::SourceLocation &source);

        ///\brief Read a source with an Elliptical Gaussian PSF that can be 
        ///enabled/disabled according to a column in the input file
        SourceInput &operator>>(IO::OutputSDKSource &source);

        ///\brief Is the state of the source reader bad (i.e. last operation 
        ///failed or we reached the end of input
        bool operator!() const {return !bool(*this);}

        ///Did the last operation complete successfully.
        bool fail() const {return !__good;}

        ///\brief Is the state of the source reader good (i.e. last operation 
        ///succeeded and no EOF yet.
        operator bool() const {return __good && !__eof;}

        ///Did we reach the end of input
        bool eof() const {return __eof;}
    }; //End SourceInput class.

} //End SubPixPhot namespace.

#endif
