#include "IOSources.h"

namespace FitPSF {
    typedef std::list<std::string> ColNameList;
    typedef std::vector< std::list<double> > ColumnData;

    void IOSources::read_filenames(std::istream &input_stream)
    {
        if(!input_stream.good())
            throw Error::IO(
                "Requested sources from a stream in non-good state"
            );

        int error_position = input_stream.tellg();

        input_stream >> __fits_fname >> __output_fname;
        if(__output_fname[__output_fname.size() - 1] == ']') {
            __source_assignment_fname = "]";
            __output_fname.resize(__output_fname.size() - 1);
        } else input_stream >> __source_assignment_fname;
        if(
            !input_stream.good()
            ||
            __fits_fname[0] != '['
            ||
            __source_assignment_fname[
                __source_assignment_fname.size() - 1
            ] != ']'
        ) {
            std::ostringstream message;
            message << "FitPSF input section (at character "
                    << error_position
                    << ") does not begin with "
                    << "[<FITS filename> <output filename> "
                    << "<source assignment filename>]";
            throw Error::IO(message.str());
        }
        __fits_fname.erase(0, 1);
        __source_assignment_fname.resize(
            __source_assignment_fname.size() - 1
        );
    }

    void IOSources::read_column_data(std::istream &input_stream,
                                     const ColNameList &column_names,
                                     ColumnData &column_data)
    {
        while(input_stream.good()) {
            while(std::isspace(input_stream.peek())) {
                if(input_stream.eof()) {
                    __last = true;
                    return;
                }
                input_stream.get();
            }
            if(input_stream.eof() || input_stream.peek() == '[') {
                __last = input_stream.eof();
                return;
            }

            std::string line;
            std::getline(input_stream, line);
            std::istringstream line_stream(line);

            int error_position = input_stream.tellg();
            std::string source_id;
            double x = Core::NaN, y = Core::NaN, value = Core::NaN;

            ColumnData::iterator destination=column_data.begin();
            for(
                ColNameList::const_iterator col_i = column_names.begin();
                col_i!=column_names.end();
                ++col_i
            ) {
                if(!line_stream.good()) {
                    std::ostringstream message;
                    message << "FitPSF input (at character "
                            << error_position
                            << ") contains a malformatted source line";
                    throw Error::IO(message.str());
                }

                if(*col_i == "ID") line_stream >> source_id;
                else {
                    line_stream >> value;
                    if(*col_i == "x") x = value;
                    else if(*col_i == "y") y = value;
                    destination->push_back(value);
                    ++destination;
                }
            }
            __locations.push_back(
                Core::SourceLocation(Core::SourceID(source_id), x, y)
            );
        }
        __last = input_stream.eof();
    }

    void fill_valarray_from_list(const std::list<double> &from,
                                 std::valarray<double> &to)
    {
        to.resize(from.size());
        size_t to_i = 0;
        for(
            std::list<double>::const_iterator from_i = from.begin();
            from_i != from.end();
            ++from_i
        ) to[to_i++] = *from_i;
    }

    IOSources::IOSources(std::istream &input_stream,
                         const ColNameList &column_names)
    {
        read_filenames(input_stream);
        std::string line;
        ColumnData column_data(column_names.size() - 1);
        read_column_data(input_stream, column_names, column_data);

        ColumnData::const_iterator column_values = column_data.begin();
        for(
            ColNameList::const_iterator col_i = column_names.begin();
            col_i!=column_names.end();
            ++col_i
        ) {
            if(*col_i != "ID") {
                __columns.push_back(
                        PSF::MapVariableType(*col_i, std::valarray<double>())
                );
                fill_valarray_from_list(*column_values,
                                        __columns.back().second);
                ++column_values;
            }
        }
    }
}
