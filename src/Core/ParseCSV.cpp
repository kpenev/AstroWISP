/**\file
 *
 * \brief The definitions of the function in ParseCSV.h
 *
 * \ingroup IO
 */

#include "ParseCSV.h"

namespace Core {

    ///\brief Extends a list with a number of copies of its last elements.
    ///
    ///Because the output contains columns of quantities that are aperture
    ///dependent (as well as some that are not), it is necessary to create
    ///multiple copies of the column identifiers that need to have one
    ///column per aperture.
    template<class ITER_TYPE, class VAL_TYPE>
    void replicate_last(
            ///All elements after this one (including this one) are
            ///replicated.
            ITER_TYPE &start_replicate,

            ///The number of copies to make.
            unsigned copies,

            ///The list this operation is performed on.
            std::list<VAL_TYPE> &target)
    {
        if(copies==0 || start_replicate==target.end()) return;
        target.push_back(*start_replicate);
        ITER_TYPE stop_replicate=target.end(), skip1=start_replicate;
        --stop_replicate;
        ++skip1;
        target.insert(target.end(), skip1, stop_replicate);
        for(;copies>1; --copies)
            target.insert(target.end(), start_replicate, stop_replicate);
    }



    RealList parse_real_list(const std::string &csv,
                             const std::string &optname,
                             unsigned min_count,
                             unsigned max_count)
    {
        RealList result;
        parse_csv_list(csv, optname, min_count, max_count, result);
        return result;
    }

    std::list<int> parse_int_list(const std::string &csv,
                                  const std::string &optname,
                                  unsigned min_count,
                                  unsigned max_count)
    {
        std::list<int> result;
        parse_csv_list(csv, optname, min_count, max_count, result);
        return result;
    }

    ColumnList parse_column_list(const std::string &csv,
                                       unsigned num_apertures,
                                       const std::string &optname,
                                       bool allow_unknown)
    {
        std::istringstream csv_stream(csv);
        ColumnList result;
        std::string colname;
        std::list<Phot::Columns>::const_iterator to_replicate=result.end();
        bool per_ap=false, old_per_ap=false;
        while(csv_stream) {
            getline(csv_stream, colname, ',');
            if(colname == "id") {
                result.push_back(Phot::id);
                per_ap = false;
            } else if(colname == "x") {
                result.push_back(Phot::x);
                per_ap = false;
            } else if(colname == "y") {
                result.push_back(Phot::y);
                per_ap = false;
            } else if(colname == "S") {
                result.push_back(Phot::S);
                per_ap = false;
            } else if(colname == "D") {
                result.push_back(Phot::D);
                per_ap = false;
            } else if(colname == "K") {
                result.push_back(Phot::K);
                per_ap = false;
            } else if(colname == "A" || colname == "amp") {
                result.push_back(Phot::A);
                per_ap = false;
            } else if(colname == "bg") {
                result.push_back(Phot::bg);
                per_ap = false;
            } else if(colname == "bg_err") {
                result.push_back(Phot::bg_err);
                per_ap = false;
            } else if(colname == "flux") {
                result.push_back(Phot::flux);
                per_ap = true;
            } else if(colname == "flux_err") {
                result.push_back(Phot::flux_err);
                per_ap = true;
            } else if(colname == "mag") {
                result.push_back(Phot::mag);
                per_ap = true;
            } else if(colname == "mag_err") {
                result.push_back(Phot::mag_err);
                per_ap = true;
            } else if(colname == "flag") {
                result.push_back(Phot::flag);
                per_ap = true;
            } else if(colname == "enabled") {
                result.push_back(Phot::enabled);
                per_ap = false;
            } else if(colname == "chi2") {
                result.push_back(Phot::chi2);
                per_ap = false;
            } else if(colname == "sn") {
                result.push_back(Phot::sn);
                per_ap = false;
            } else if(colname == "npix") {
                result.push_back(Phot::npix);
                per_ap = false;
            } else if (colname == "nbgpix") {
                result.push_back( Phot::nbgpix );
                per_ap = false;
#ifdef DEBUG
            } else if (colname == "time") {
                result.push_back(Phot::time);
                per_ap = false;
#endif
            } else if(allow_unknown) {
                result.push_back(Phot::unknown);
                per_ap = false;
            } else if(colname != "") {
                std::ostringstream msg;
                msg << "Unrecognized column '" << colname << "'"
                    <<  " in " << optname;
                throw Error::CommandLine(msg.str());
            }
            if (old_per_ap && !per_ap && num_apertures > 1) {
                replicate_last(to_replicate, num_apertures - 1, result);
            }
            if(per_ap && !old_per_ap) {
                to_replicate = result.end();
                --to_replicate;
            }
            if(csv_stream.eof()) {
                if (per_ap && num_apertures > 1) {
                    replicate_last(to_replicate, num_apertures - 1, result);
                }
                return result;
            }
            old_per_ap = per_ap;
        }
        std::ostringstream error_msg;
        error_msg << "Malformatted " << optname << " option: " << csv
                  << "expected comma separated list of column names.";
        throw Error::CommandLine(error_msg.str());
    }

    StringList parse_string_list(const std::string &csv,
                                 const std::string &optname,
                                 unsigned min_count,
                                 unsigned max_count)
    {
        std::istringstream csv_stream(csv);
        StringList result;
        while(csv_stream) {
            result.push_back("");
            getline(csv_stream, result.back(), ',');
            if(
                csv_stream.eof()
                &&
                result.size()>=min_count
                &&
                result.size()<=max_count
            )
                return result;
        }
        std::ostringstream error_msg;
        error_msg << "Malformatted " << optname << " option: " << csv
                  << "expected comma separated list of "
                  << "at least " << min_count << " "
                  << "and "
                  << "at most " << max_count << " "
                  << "real values.";
        throw Error::CommandLine(error_msg.str());
    }

}
