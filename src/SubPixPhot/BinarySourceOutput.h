/**\file
 *
 * \brief Declare a class for adding sources to a binary photometry file.
 *
 * \ingroup SubPixPhot
 */

#ifndef __BINARY_SOURCE_OUTPUT_H
#define __BINARY_SOURCE_OUTPUT_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Core/PhotColumns.h"
#include "../IO/OutputSDKSource.h"
#include "../IO/binostream.h"

namespace SubPixPhot {

    ///\brief A class that behaves as output stream for sources resulting in 
    ///binary photometry file.
    class LIB_LOCAL BinarySourceOutput {
    private:
        ///The stream where output will be written
        IO::binostream *__outstream;

        ///\brief The columns to output in the order in which they should be 
        ///output.
        ///
        ///Aperture specific columns may appear as many times as 
        ///there are apertures.
        const std::list<Phot::Columns> *__columns;


        ///\brief Indexed by the value of Phot::Columns corresponding to the 
        ///given column, this array specifies which columns are necessary on 
        ///output
        std::valarray<bool> __required;

        ///\brief The presicion with which to output each column (indexed 
        ///by the  value of Phot::Columns). Entries corresponding to
        ///integer columns are ignored.
        std::valarray<double> __precision;

        ///How many apertures are requried by the output columns
        unsigned __num_apertures;

        ///What magnitude corresponds to a flux of 1ADU
        double __mag_1ADU;

        ///\brief Fills in the __required array according to the content of
        ///__columns.
        void set_required_columns();

        ///\brief Returns the format string that encodes what non-photometry 
        ///columns are in the file.
        std::string aperture_indep_format();

        ///\brief Retruns the format string that encodes what photometry 
        ///columns are in the file.
        std::string aperture_dep_format();
    public:
        ///Use the given stream for output.
        BinarySourceOutput(
            ///The stream to write the sources to.
            IO::binostream &outstream,

            ///The magnitude that corresponds to a flux of 1ADU.
            double mag_1ADU,

            ///The columns to include in the output.
            const std::list<Phot::Columns> &output_columns =
                std::list<Phot::Columns>()
        ); 

        ///Sets the columns (including order) which to output.
        void set_columns(
            ///The new columns to include in the output.
            const std::list<Phot::Columns> &output_columns
        );

        ///\brief Outputs the given list of sources to the given stream in 
        ///packed binary format.
        void operator()(
            ///The sources which were just photometered to output.
            const std::list<IO::OutputSDKSource> &sources
        );

        ///The precision with which the given column will be output.
        double precision(
            ///The column to get the precision of.
            Phot::Columns column
        ) const
        {return __precision[column];}

        ///The precision with which the given column will be output.
        double &precision(
            ///The column to get the precision of.
            Phot::Columns column
        )
        {return __precision[column];}
    }; //End BinarySourceOutput class.

} //End SubPixPhot namespace.

#endif
